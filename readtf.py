from Bio import SeqIO
import gzip
import sys
import os
import pymysql

hostname = 'localhost'
username = 'root'
password = ''
database = 'bimm185'

def main():
	myConnection = pymysql.connect(host=hostname, user=username, passwd=password, db=database, local_infile=True, autocommit=True)

	with open("network_tf_gene.txt",'r') as file:
		for line in file:
			if "#" in line:
				continue
			else:
				line = line.strip().split('\t')
				print(line)






if __name__ == '__main__':
	main()