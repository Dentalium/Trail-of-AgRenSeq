#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO, Seq
import gzip
import io
import sys
from Phenotype_GLM import Phenotype
import pandas as pd
from sklearn.decomposition import PCA
import statsmodels.formula.api as smf
import statsmodels.api as sm # debug
import numpy as np
import math
import time
from bitarray import bitarray 


"""
Correlation pre-filtering of k-mers followed by 
(i) regression analysis to control for population structure using PCA, and 
(ii) projection of k-mers onto contigs of an assembly.

Negative log of p-values obtained after regression analysis are stored as the measures of association.

@authors: kumar gaurav & sanu arora
"""
class KmerProjection(object):
    
    """
    @param phenotype
    @param kmerSize
    @param assemblyFile: fasta file containing the denovo assembly of the accession onto which k-mers are projected and their association scores stored. 
    @param nlrList: tsv file with contig ids in the first column. Only these contigs are retained from the denovo assembly.
    @param presenceMatrix: 
                    Gzipped text file. Unzipped not supported. 
                    First line is assumed to start with a # followed by a comma separated list of accession names. 
                    An entry in the matrix is a kmer, then a tab, then a string of 0 and 1 for absence and presence of that kmer in accessions. Order is according to first line of the matrix. Trailing 0s are ommitted.
    @param snpFile: tab separated file containing SNP markers presence/absence matrix
    @param pcaDimensions: number of significant PCA dimensions retained for regression analysis
    @param cor_threshold: correlation threshold for pre-filtering
    """
    def __init__(self, phenotype, assemblyFile, nlrList, presenceMatrix, snpFile, pcaDimensions, cor_threshold):

        self.associationMatrix = {}
        self.phenotype = phenotype
        self.cor_threshold = cor_threshold
        
        self.assemblyFile = assemblyFile
        self.nlrList = nlrList
        self.presenceMatrix = presenceMatrix
        self.snpFile = snpFile
 
        self.kmerSize = self.setKmerSize()
        
        self.pcaDF = self.computePCA(pcaDimensions)

        # predict phenotype scores using only PCA dimensions
        self.nullDF = self.pcaDF.copy()
        self.nullDF['bias'] = np.ones(self.nullDF.shape[0]).reshape(-1,1)
        # statsmodels 0.11.1不支持smf.OLS()的写法
        #self.null_results = smf.OLS(self.phenotype.phenoScores_series, self.nullDF).fit()
        self.null_results = sm.OLS(self.phenotype.phenoScores_series, self.nullDF).fit()


    """
    @param n_dimensions
    """    
    def computePCA(self, n_dimensions):

        snpDF = pd.DataFrame.from_csv(self.snpFile, sep = "\t")

        ## debug
        #print(self.phenotype.phenoScores_series.index)
        #print(snpDF.index)

        # retain only those accessions in SNP markers matrix for which phenotype scores are stored
        snpDF_reduced = snpDF.loc[self.phenotype.phenoScores_series.index]

        print("Computing PCA with " + str(n_dimensions) + " components.")
        pca = PCA(n_components = n_dimensions).fit_transform(snpDF_reduced)
        pcaDF= pd.DataFrame(pca, index = snpDF_reduced.index)
        pcaDF.columns = ['PCA_dim_' + str(col+1) for col in pcaDF.columns]

        return pcaDF


    """
    calculate correlation.     
    @param presenceString: Ordered string of 0 and 1 representing presence and absence of k-mer in the accessions. 
    @param accessions_dict: mapping of accessions to their index in the presnceString. 
    @params pheno_sum, rho_pheno: phenotype parameters used in correlation calculation.
    @param n_accessions: number of accessions in accessions_dict. explicitly passed to reduce computation time.
    """
    def getCorrelation(self, presenceString, accessions_dict, pheno_sum, rho_pheno, n_accessions):
        
        presence_bitarray = bitarray(presenceString)
        
        presence_sum = 0.0
        presence_sqsum = 0.0
        presence_pheno_dot_product = 0.0

        # debug
        '''
        print("length of accession", len(self.phenotype.phenoScores_dict))
        print(list(self.phenotype.phenoScores_dict.items())[:5])
        print("accessions_dict", len(accessions_dict))
        print("accessions_dict 内容示例:", list(accessions_dict.items())[:5])
        '''

        """
        self.phenotype.phenoScores_dict: 一个字典, 储存了表型(抗性得分)
        accessions_dict: 一个字典, 储存了样品在kmer矩阵中的位置(索引, 从0开始), 长度为151
        """

        for accession in self.phenotype.phenoScores_dict:
            index = accessions_dict[accession]
            # debug
            '''
            print("presence_bitarray 长度:", len(presence_bitarray))
            print("当前 index:", index)
            '''
            # 加上try语句避免由于kmer矩阵样品数大于表型文件样品数导致的索引越界
            try:
                if presence_bitarray[index]:
                    presence_sum += 1
                    presence_sqsum += 1
                    presence_pheno_dot_product += self.phenotype.phenoScores_dict[accession]
            except:
                print(accession, "excluded")

        rho_presence = math.sqrt(n_accessions*presence_sqsum - presence_sum*presence_sum)

        # if rho_presence is 0, k-mer is either present in all the accessions or absent in all the accessions, therefore not informative for association analysis.
        # to avoid division by 0 in correlation calculation, k-mer is given an association score of -1.0, ensuring that it is pre-filtered.
        if rho_presence == 0:
            return -1.0 

        correlationScore = (n_accessions*presence_pheno_dot_product - presence_sum*pheno_sum)/(rho_presence*rho_pheno)    

        return correlationScore

        
    """
    regression analysis.     
    @param presenceString: Ordered string of 0 and 1 representing presence and absence of k-mer in the accessions of header. 
    @param header 
    """    
    def getGLMpval(self, presenceString, header):
       
        # debug
        print(presenceString, "len:", len(presenceString))
        print(header, "len:", len(header))

        # 补0
        if len(presenceString)<len(header):
            print("filling tailing 0")
            presenceString+="0"*(len(header)-len(presenceString))

        presence_array = np.array(list(presenceString)).astype(int)
        presence_series = pd.Series(presence_array, index = header)

        pca_presenceDF = self.nullDF.copy()
        pca_presenceDF['presence'] = presence_series.loc[self.phenotype.phenoScores_series.index]
        
        # predict phenotype scores based on presence/absence of k-mer, with PCA dimensions as covariates        
        # 同上，改为sm.OLS
        #results = smf.OLS(self.phenotype.phenoScores_series, pca_presenceDF).fit()    
        results = sm.OLS(self.phenotype.phenoScores_series, pca_presenceDF).fit()

        # pvalue for nested models
        pval = results.compare_lr_test(self.null_results)[1]

        return -np.log(pval)
    

    """
    read through the file containing the presence/absense strings of kmers, pre-filter k-mers based on correlation threshold, calculate association scores
    of pre-filtered k-mers based on regression analysis.
    """    
    def readMatrix_GLM(self):
        
        print("Reading presence/absense matrix " + self.presenceMatrix)
        gz = gzip.open(self.presenceMatrix, 'r')
        f = io.BufferedReader(gz)
        
        accessions_dict = {}
        header = next(f).decode("utf-8") 
        header = [x.strip() for x in header.split(',')]     
        header[0] = header[0][1:]
        
        for i in range(len(header)):
            accessions_dict[header[i]] = i
        
        phenoScores_dict = self.phenotype.phenoScores_dict
        print("Correlation pre-filtering with a threshold of " + str(self.cor_threshold))

        # calculate phenotype parameters used in correlation calculation. reduces computation time (not by much).
        pheno_sum = 0.0
        pheno_sqsum = 0.0
    
        for accession in phenoScores_dict:
            pheno_sum += phenoScores_dict[accession]
            pheno_sqsum += phenoScores_dict[accession]*phenoScores_dict[accession]
            
        n_accessions = len(phenoScores_dict)
        rho_pheno = math.sqrt(n_accessions*pheno_sqsum - pheno_sum*pheno_sum)

        # rho_pheno is 0 if and only if phenotype scores are all equal. In that case, no association mapping is possible, therefore exit the program. 
        if rho_pheno == 0:
            print("Phenotype scores used in the correlation calculation cannot all be same.")
            sys.exit()    

        # initialize parameters to track progress 
        progress = 0
        t_init = time.time()

        for inputline in f:
            # debug
            """print(inputline, len(inputline.decode()))"""
            # track progress
            progress += 1
            if progress%1000000 == 0:
                print (str(progress)+" k-mers parsed in time: " +str(time.time() - t_init))

            split = inputline.decode("utf-8").split()
            kmerF = split[0].upper()
            kmerR = str(Seq.Seq(kmerF).reverse_complement())
            bool_kmerF = kmerF in self.associationMatrix
            bool_kmerR = kmerR in self.associationMatrix

            if bool_kmerF:
                correlation = self.getCorrelation(split[1], accessions_dict, pheno_sum, rho_pheno, n_accessions)
                if correlation > 1.0 or correlation < -1.0:
                    print(correlation)

                if correlation < self.cor_threshold:
                    # assign a negative p-value to pre-filter
                    pvalF = -1.0    
                else:
                    #print(progress)
                    pvalF = self.getGLMpval(split[1], header)

                self.associationMatrix[kmerF] = pvalF  

                if bool_kmerR:
                    self.associationMatrix[kmerR] = pvalF 

            else:
                if bool_kmerR:
                    correlation = self.getCorrelation(split[1], accessions_dict, pheno_sum, rho_pheno, n_accessions)
                    if correlation < self.cor_threshold:
                        pvalR = -1.0    
                    else:
                        pvalR = self.getGLMpval(split[1], header)

                    self.associationMatrix[kmerR] = pvalR  
                                        
        gz.close()
    

    """
    computes the k-mer length from the first entry in the presence absense matrix.
    """
    def setKmerSize(self):

        gz = gzip.open(self.presenceMatrix, 'r')
        f = io.BufferedReader(gz)
        next(f)
        inputline = f.readline()
        gz.close()
        return len(inputline.split()[0])
        
        
    """
    k-mers extracted from a denovo assembly, and only those retained which also appear in the list of NLR contigs 
    nlrList is assumed to be a tab separated table, containing the contig ids in the first column.
    """
    def readAssembly_GLM(self):
        
        print("Reading assembly " + self.assemblyFile + " using kmer size of " + str(self.kmerSize))
        
        nlrContigs = set()
        with open(self.nlrList) as f:
            next(f)
            for l in f.readlines():
                nlrContigs.add(l.split()[0])

        print("number of NLRs: "+str(len(nlrContigs)))
        
        fasta_sequences = SeqIO.parse(open(self.assemblyFile),'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            if name not in nlrContigs:
                continue
            for i in range(len(sequence) - self.kmerSize + 1):
                kmer =  sequence[i:i+self.kmerSize].upper()
                # possible to encode kmers as bitarrays, but no practically significant difference in resource requirements.
                self.associationMatrix[kmer] = -1.0

        print("...finished. Recorded " + str(len(self.associationMatrix)) + " kmers.")
    

    """
    k-mers mapped to NLR contigs of denovo assembly, and their association scores recorded, if positive.     
    @param outputFile: tsv file
                First column - contig ID. 
                Second column - a sequence of consecutive natural numbers along the first column. Used as x-axis for plotting association scores.
                Third column - unique associations scores of k-mers mapping to that contig.
                Fourth column - number of kmers having that score.
    """
    def writeAssociationScore_GLM(self, outputFile):
		
        print("Writing association output")
        
        nlrContigs = set()
        with open(self.nlrList) as f:
            next(f)
            for l in f.readlines():
                nlrContigs.add(l.split()[0])
                
        fasta_sequences = SeqIO.parse(open(self.assemblyFile),'fasta')
        out = open(outputFile, 'w')
        contigCount = 0 
        
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            if name not in nlrContigs:
                continue
            contigCount += 1
            h = {}
            for i in range(len(sequence) - self.kmerSize + 1):                
                kmer = sequence[i:i+self.kmerSize].upper()
                associationScore  = self.associationMatrix[kmer] 
                if associationScore > 0.0:
                    num = 0
                    if associationScore in h:
                        num = h[associationScore]
                    num += 1
                    h[associationScore] = num
                                 
            for associationScore, num in h.items():
                out.write(name + "\t" + str(contigCount) + "\t" + str(associationScore) + "\t" + str(num) + "\n")

        out.close()

    

        

    
        
