#This code calculates the rho and distance by:
#1. querying gene pairs from the database
#2. reading the expression data from colombos file
#3. need to uncomment the operon reading part to store the operon info into the database


import gzip
import sys
import os
import pymysql
import scipy
from scipy import stats
import matplotlib.pyplot as plt
import numpy
from numpy import linspace
from scipy.stats.kde import gaussian_kde


hostname = 'localhost'
username = 'root'
password = ''
database = 'bimm185'


def query_gene(conn, gene):
	cur=conn.cursor()
	sql_statement = ("SELECT gene_id, name FROM GENES WHERE name='{gene}' and genome_id = 1".format(gene=gene))
	cur.execute(sql_statement)
	result = cur.fetchone()
	return result

def query_gene_locau_tag(conn, locus_tag):
	cur = conn.cursor()
	sql_statement = ("SELECT gene_id, name FROM genes WHERE locus_tag = '{locus_tag}' and genome_id = 1".format(locus_tag=locus_tag))
	cur.execute(sql_statement)

	result = cur.fetchone()

	if result is None:
		return None
	else:
		return result

def query_directons(conn):
	cur = conn.cursor()
	sql_statement = ("select genes.gene_id, genes.name, strand, left_position, right_position from genes inner join(select gene_id, min(left_position) as left_position, max(right_position) as right_position from exons group by gene_id) position on position.gene_id = genes.gene_id where genes.genome_id = 1 order by left_position;")
	cur.execute(sql_statement)
	result = cur.fetchall()
	#print(result)
	return result

def query_synonyms(conn, gene):
	cur=conn.cursor()
	sql_statement = ("SELECT genes.gene_id, name FROM gene_synonyms inner join genes on gene_synonyms.gene_id = genes.gene_id WHERE synonym='{gene}' and genes.genome_id = 1".format(gene=gene))
	cur.execute(sql_statement)
	result = cur.fetchone()
	return result

def import_tf_gene(conn):
	cur = conn.cursor()
	sql_statement = ("LOAD DATA LOCAL INFILE 'tf_gene.txt' INTO TABLE TF_gene;")
	cur.execute(sql_statement)


def create_table(conn):
	cur = conn.cursor()

	sql_statement = ("drop table if exists TF_gene;"
	"CREATE TABLE TF_gene (\n"
	"TF VARCHAR(256) NOT NULL,\n"
	"gene_id INT (10) UNSIGNED NOT NULL,\n"
	"gene_name VARCHAR(256) NOT NULL,\n"
	"effect VARCHAR(256) NOT NULL\n"
	")ENGINE=InnoDB;")
	#print(sql_statement)
	cur.execute(sql_statement)
	cur.close()


def create_operon_table(conn):
	cur = conn.cursor()
	sql_statement = ("drop table if exists operons;"
	"CREATE TABLE operons (\n"
	"operon_id INT (10) NOT NULL,\n"
	"gene_id INT (10) UNSIGNED NOT NULL\n"
	")ENGINE=InnoDB;")
	cur.execute(sql_statement)
	cur.close()

def import_operon(conn):
	cur = conn.cursor()
	sql_statement = ("LOAD DATA LOCAL INFILE 'operons.txt' INTO TABLE operons;")
	#print(sql_statement)
	cur.execute(sql_statement)
	#print(sql_statement)
	cur.close()

def create_colombos(conn):

	cur = conn.cursor()
	sql_statement = ("drop table if exists colombos;"
	"CREATE TABLE colombos (\n"
	"gene_id INT (10) NOT NULL, \n"
	"expression double precision, \n"
	"c VARCHAR(256) NOT NULL \n"
	")ENGINE=InnoDB;")
	#print(sql_statement)
	cur.execute(sql_statement)
	cur.close()

def import_colombos(conn):
	cur = conn.cursor()
	sql_statement = ("LOAD DATA LOCAL INFILE 'colombos.txt' INTO TABLE colombos;")
	#print(sql_statement)
	cur.execute(sql_statement)
	#print(sql_statement)
	cur.close()

def alter_colombos(conn, attrs):
	cur = conn.cursor()
	for a in attrs:
		sql_statement = ("ALTER TABLE colombos ADD c{a} DOUBLE;".format(a=a))
		cur.execute(sql_statement)
	cur.close()

def query_gene_pairs(conn):
	cur = conn.cursor()
	sql_statement=("select t1.gene_id, t2.gene_id from (select * from TF_gene where TF not in (select TF from (select count(*) as count, TF from TF_gene group by TF) t where count > 130)) t1, TF_gene t2 where t1.gene_id < t2.gene_id and t1.TF = t2.TF and t1.effect = t2.effect group by t1.gene_id, t2.gene_id;")
	cur.execute(sql_statement)
	result = cur.fetchall()
	cur.close()
	return result

def query_expressions(conn, g1, g2):
	cur = conn.cursor()
	sql_statement=("select c1.expression, c2.expression from (select * from colombos where gene_id = {g1}) c1 inner join (select * from colombos where gene_id = {g2}) c2 on c1.c = c2.c ".format(g1=g1, g2=g2))
	cur.execute(sql_statement)
	result = cur.fetchall()
	cur.close()
	return result

def create_directons_table(conn):
	cur = conn.cursor()
	sql_statement = ("drop view if exists directons;create view directons as (select genes.gene_id, genes.name, strand, left_position, right_position from genes inner join(select gene_id, min(left_position) as left_position, max(right_position) as right_position from exons group by gene_id) position on position.gene_id = genes.gene_id where genes.genome_id = 1 order by left_position);")
	cur.execute(sql_statement)
	result = cur.fetchall()
	#print(result)
	return result

def query_gene_position(conn, g1, g2):
	cur = conn.cursor()
	sql_statement=("select left_position, right_position from directons where gene_id = {g1} or gene_id = {g2} order by left_position".format(g1=g1, g2=g2))
	cur.execute(sql_statement)
	result = cur.fetchall()
	return result

def query_operon_pairs(conn):
	cur = conn.cursor()
	sql_statement = ("select o1.gene_id ,o2.gene_id from operons o1, operons o2 where o1.operon_id = o2.operon_id and o1.gene_id < o2.gene_id")
	cur.execute(sql_statement)
	result = cur.fetchall()
	return result


def read_gene_products():
	gene_products = {}
	with open("GeneProductSet.txt",'r') as file:
		for line in file:
			if '#' in line:
				continue
			else:
				line = line.strip().split('\t')
				if len(line) < 3:
					continue
				#print(line_count, "length of line: " , len(line))
				gene_products[line[1]] = line[2]
	return gene_products

def read_tf_gene(myConnection, gene_products):
	create_table(myConnection)
	tf_gene_file = open("tf_gene.txt",'w')

	with open("network_tf_gene.txt",'r') as file:
		for line in file:
			if "#" in line:
				continue
			else:
				line = line.strip().split('\t')
				if line[-1] == "Strong"  or len(line[-2].split(',')) > 1:
					result = query_gene(myConnection,line[1])
					if result is None:
						if line[1] in gene_products:
							locus_tag = gene_products[line[1]]
							result = query_gene_locau_tag(myConnection, locus_tag)
						if result is None:
							result = query_synonyms(myConnection, line[1])

					if result is not None:
						print(result[0],result[1], sep="\t")
						tf_gene_file.write(line[0]+"\t"+str(result[0])+"\t"+result[1]+'\t'+ line[2]+"\n")

	import_tf_gene(myConnection)
	tf_gene_file.close()

def read_operons(myConnection):
	create_operon_table(myConnection)
	operon_count = 0
	operon_file = open("operons.txt",'w')
	directons = query_directons(myConnection)
	h0_files = open("h0.txt",'w')
	h0_pairs = []
	with open("OperonSet.txt",'r') as file:
		for line in file:
			if "#" in line:
				continue
			line = line.strip().split('\t')
			if len(line) != 8:
				continue
			gene_name = line[5]
			confidence = line[7]
			strand = line[3]
			if confidence == 'Strong' or confidence == 'Confirmed':
				gene_name = gene_name.split(',')
				if len(gene_name) < 2:
					continue
				operon_count += 1
				for g in gene_name:
					if "<" in g:
						continue
					#print(g)
					result = query_gene(myConnection,g)
					if result is None:
						if line[1] in gene_products:
							locus_tag = gene_products[g]
							result = query_gene_locau_tag(myConnection, locus_tag)
						if result is None:
							result = query_synonyms(myConnection, g)

					if result is not None:
						operon_file.write(str(operon_count) + '\t' + str(result[0]) + '\n')
				if strand == 'forward':
					#first gene is the leftmost
					for i in range(len(genes)-1):
						operons.add((genes[i],genes[i+1]))

					if genes[0] in sorted_genes:
						gene_index = sorted_genes.index(genes[0])
						
						if gene_index != 0:

							prev = directons[gene_index - 1]

							#print(prev)
							if prev[0] not in border_dic:
								if prev[2] == '+':
									distance = directons[gene_index][3] - prev[4]
									borders.add((prev[0],genes[0]))
									h0.append(distance)
									border_dic[prev[0]] = genes[0]
									border_dic[genes[0]] = prev[0]
									print(prev[0],genes[0])
								else:
									print(prev[0],genes[0],'not on the same strand')
							elif border_dic[prev[0]] != genes[0]:
								if prev[2] == '+':
									distance = directons[gene_index][3] - prev[4]
									borders.add((prev[0],genes[0]))
									h0.append(distance)
									border_dic[prev[0]] = genes[0]
									border_dic[genes[0]] = prev[0]
									print(prev[0],genes[0])
								else:
									print(prev[0],genes[0],'not on the same strand')

								#print(directons[gene_index],prev,distance)

					if genes[-1] in sorted_genes:
						gene_index = sorted_genes.index(genes[-1])
						if gene_index != len(sorted_genes) - 1:
							next_gene = directons[gene_index + 1]
							borders.add((genes[-1],next_gene[0]))
							if next_gene[0] not in border_dic:
								if next_gene[2] == "+":
									distance = next_gene[3] - directons[gene_index][4]
									h0.append(distance)
									border_dic[next_gene[0]] = genes[-1]
									border_dic[genes[-1]] = next_gene[0]
									print(next_gene[0],genes[-1])
								else:
									print(next_gene[0],genes[-1], 'not on the same strand')
							elif border_dic[next_gene[0]] != genes[-1]:
								if next_gene[2] == "+":
									distance = next_gene[3] - directons[gene_index][4]
									h0.append(distance)
									border_dic[next_gene[0]] = genes[-1]
									border_dic[genes[-1]] = next_gene[0]
									print(next_gene[0],genes[-1])
								else:
									print(next_gene[0],genes[-1],'not on the same strand')
								#print(directons[gene_index],next_gene,distance)


				else:
					for i in range(len(genes)-1,-1,-1):
						operons.add((genes[i],genes[i-1]))

					if genes[0] in sorted_genes:
						gene_index = sorted_genes.index(genes[0])
						if gene_index != len(sorted_genes) - 1:

							next_gene = directons[gene_index + 1]
							borders.add((genes[0],next_gene[0]))
							if next_gene[0] not in border_dic:
								if next_gene[2] == "-":
									distance = next_gene[3] - directons[gene_index][4]
									h0.append(distance)
									border_dic[next_gene[0]] = genes[-1]
									border_dic[genes[-1]] = next_gene[0]
									print(next_gene[0],genes[-1])
								else:
									print(next_gene[0],genes[-1],'not on the same strand')
							elif border_dic[next_gene[0]] != genes[-1]:
								if next_gene[2] == "-":
									distance = next_gene[3] - directons[gene_index][4]
									h0.append(distance)
									border_dic[next_gene[0]] = genes[-1]
									border_dic[genes[-1]] = next_gene[0]
									print(next_gene[0],genes[-1])
								else:
									print(next_gene[0],genes[-1],'not on the same strand')

								#print(directons[gene_index],next_gene,distance)
						

					if genes[-1] in sorted_genes:
						gene_index = sorted_genes.index(genes[-1])
						if gene_index != 0:
							prev = directons[gene_index - 1]
							borders.add((prev[0],genes[-1]))
							#print(prev)
							if prev[0] not in border_dic:
								if prev[2] == '-':
									distance = directons[gene_index][3] - prev[4]
									h0.append(distance)
									border_dic[prev[0]] = genes[0]
									border_dic[genes[0]] = prev[0]
									print(prev[0],genes[0])
								else:
									print(prev[0],genes[0],'not on the same strand')
							elif border_dic[prev[0]] != genes[-1]:
								if prev[2] == '-':
									distance = directons[gene_index][3] - prev[4]
									h0.append(distance)
									border_dic[prev[0]] = genes[0]
									border_dic[genes[0]] = prev[0]
									print(prev[0],genes[0])
								else:
									print(prev[0],genes[0],'not on the same strand')


				
	import_operon(myConnection)

	operon_file.close()

def read_colombos(myConnection):
	create_colombos(myConnection)
	
	#colombos_file = open("colombos.txt",'w')
	colombos = {}
	with open("colombos_ecoli_exprdata_20151029.txt",'r') as file:

		for i in range(6):
			file.readline()
		header = file.readline()
		attrs = header.strip().split('\t')[3:]
		#alter_colombos(myConnection, attrs)
		for line in file:
			line = line.strip().split('\t')
			b_num = line[0]
			scores = line[3:]
			result = query_gene_locau_tag(myConnection,b_num)
			if result is not None:
				colombos[result[0]] = scores
				#for i in range(len(attrs)):
				#	if scores[i] == "NaN":
				#		continue
					#colombos_file.write(str(result[0])+'\t' + scores[i] + '\t' + attrs[i] +'\n')

	
	#import_colombos(myConnection)
	#colombos_file.close()
	return colombos

def pdf_calculate(h0, h1):
	dist_space0 = linspace(-1, 1, 1000)

	kde1 = gaussian_kde(h1, 0.5)
	#kde1.set_bandwidth(bw_method='silverman')
	#kde1.set_bandwidth(bw_method=kde1.factor * 3.0)
	#kde1.covariance_factor = lambda : 0.25
	#kde1._compute_covariance()


	density1 = kde1(dist_space0)
	max_density1 = max(density1)
	for i in range(len(density1)):
		density1[i] = density1[i]/max_density1
	plt.plot(dist_space0, kde1(dist_space0)/max(kde1(dist_space0)))
	

	kde0 = gaussian_kde(h0,0.5)
	#kde0.set_bandwidth(bw_method=kde0.factor * 3.0)
	#kde0.covariance_factor = lambda : 0.25
	#kde0._compute_covariance()
	
	density0 = kde0(dist_space0)
	max_density0 = max(density0)
	dist0 = [i for i in range(min(h0),min(h0))]
	for i in range(len(density0)):
		density0[i] = density0[i]/max_density1
		#print(max_density0)
	plt.plot(dist_space0, kde0(dist_space0)/max(kde1(dist_space0)))
	plt.show()
	#print("max density0: ", max_density0)
	#print(density0)


	#print(kde0(0))
	#print(kde0(1))
	#print(kde1(0))
	#print(kde1(1))

	#plt.show()
	#return kde0, kde1



def main():
	myConnection = pymysql.connect(host=hostname, user=username, passwd=password, db=database, local_infile=True, autocommit=True)
	
	#uncomment these parts to read the data if needed
	'''
	gene_products = read_gene_products()
	read_tf_gene(myConnection)
	read_operons(myConnection)
	read_colombos(myConnection)
	'''
	#gene_products = read_gene_products()
	#read_tf_gene(myConnection, gene_products)
	colombos = read_colombos(myConnection)
	

	'''
	gene_pairs = query_gene_pairs(myConnection)
	rhos = []
	unique_rhos = []
	distances = []
	pairs = {}
	for gp in gene_pairs:
		g1 = gp[0]
		g2 = gp[1]
		e1=colombos[g1]
		e2=colombos[g2]
		#expressions = query_expressions(myConnection, g1, g2)
		#e1=[]
		#e2=[]
		positions = query_gene_position(myConnection, g1, g2)
		left = positions[1][0]
		right = positions[0][0]
		distance = left - right
		distance = min(4641652-distance, distance)
		#distances.append(distance)
		#for e in expressions:
		#	e1.append(e[0])
		#	e2.append(e[1])

		rho, pvalue = scipy.stats.spearmanr(e1, e2)
		if (g2, g1) not in pairs:
			pairs[(g1, g2)] = [rho, distance]
		else:
			pairs[(g2, g1)][0] = (pairs[(g2, g1)][0] * rho) ** 0.5


		#rhos.append(rho)
		if rho not in unique_rhos:
			unique_rhos.append(rho)

		#print(distance, rho)
		#plt.plot(distance, rho)
	#unique_rhos = sorted(unique_rhos)

	#ranks = []
	#for i in range(10000):
	#	index = unique_rhos.index(rhos[i]) + 1
	#	ranks.append(index)

	#for i in range(10000):
	#	print(distance[i], ranks[i])

	for key in pairs.keys():
		print(pairs[key][1], pairs[key][0], sep="\t")



	
	#plt.plot(distances, rhos,'ro')
	#plt.show()
	'''

	#Test with operon pairs
	
	means = {}
	gene_pairs = query_operon_pairs(myConnection)
	distances = []
	rhos = []
	for gp in gene_pairs:
		g1 = gp[0]
		g2 = gp[1]
		e1=colombos[g1]
		e2=colombos[g2]
		#expressions = query_expressions(myConnection, g1, g2)
		#e1=[]
		#e2=[]
		positions = query_gene_position(myConnection, g1, g2)
		left = positions[1][0]
		right = positions[0][0]
		distance = left - right
		distance = min(4641652-distance, distance)
		distances.append(distance)
		#for e in expressions:
		#	e1.append(e[0])
		#	e2.append(e[1])

		rho, pvalue = scipy.stats.spearmanr(e1, e2)


		rhos.append(rho)
		print(distance, rho)
		if distance in means:
			means[distance].append(rho)
		else:
			means[distance] = [rho]

	

		 #print(distance, rho)
		#plt.plot(distance, rho)
	#unique_rhos = sorted(unique_rhos)
	#plt.plot(distances, rhos, "ro")
	#plt.show()

	#calculate mean

	'''
	ds = []
	ms = []
	for key in means.keys():
		ds.append(key)
		ms.append(sum(means[key])/len(means[key]))
	'''

	'''
	means = {}
	with open("operons_rho.out", "r") as file:
		for line in file:
			line = line.strip().split(" ")
			if line[0] in means:
				means[line[0]].append(float(line[1]))
			else:
				means[line[0]] = [float(line[1])]
	ds = []
	ms = []
	for key in means.keys():
		ds.append(key)
		ms.append(sum(means[key])/len(means[key]))

	plt.plot(ds, ms, "ro")
	plt.show()
	'''


	pdf_calculate(distances, rhos)














		










if __name__ == '__main__':
	main()