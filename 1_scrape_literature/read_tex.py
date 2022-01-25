import re
import numpy as np
import os
import arxiv
import argparse
import sys

parser = argparse.ArgumentParser(description='Arxiv')
parser.add_argument('--star', '-s')
parser.add_argument('--pl_id', '-id')
args = parser.parse_args()

star = args.star
pl_id = args.pl_id

planet = star.upper()+ '-' + pl_id 
planet_name_variant_1 =  star.upper() + '-00' + pl_id
planet_name_variant_2 = star + '-' + pl_id 
planet_name_variant_3 = star + '-00' + pl_id

planet_a = star.upper()+ '-' + pl_id + 'b'
planet_name_variant_1a =  star.upper() + '-00' + pl_id + 'b'
planet_name_variant_2a = star + '-' + pl_id + 'b'
planet_name_variant_3a = star + '-00' + pl_id + 'b'


planet_name_variant_4 = star + '-' + pl_id + 'A'
planet_name_variant_4a = star + '-' + pl_id + 'A' +'b'
planet_name_variant_4b = star.upper() + '-' + pl_id + 'A'
planet_name_variant_4c = star.upper() + '-' + pl_id + 'A' +'b'

direct = os.path.dirname(os.getcwd())


IDs = [] # arXiv IDs with mid-transtit times 
count = 0
total_count = 0
for arxiv_id in os.listdir(direct + f'/arxiv/untar_dir/{planet}'):
	tex_files = 0
	for filename in os.listdir(direct + f'/arxiv/untar_dir/{planet}/{arxiv_id}'):
		if filename.endswith(".tex"):
			#print('filename: ', filename)
			data = open(direct + f'/arxiv/untar_dir/{planet}/{arxiv_id}/' + filename, errors='ignore').read()
			cond = data.count(planet) > 1 or data.count(planet_name_variant_1) > 1  or data.count(planet_name_variant_2) > 1  or data.count(planet_name_variant_3) > 1  or data.count(planet_a) > 1  or data.count(planet_name_variant_1a) > 1 or data.count(planet_name_variant_2a) > 1 or data.count(planet_name_variant_3a) > 1 or data.count(planet_name_variant_4) > 1 or data.count(planet_name_variant_4a) > 1
			cond1 = data.find(planet) != -1 or data.find(planet_name_variant_1) != -1  or data.find(planet_name_variant_2) != -1 or data.find(planet_name_variant_3) != -1 or data.find(planet_a) != -1  or data.find(planet_name_variant_1a) != -1 or data.find(planet_name_variant_2a) != -1 or data.find(planet_name_variant_3a) != -1 or data.find(planet_name_variant_4) != -1 or data.find(planet_name_variant_4a) != -1
	
			if cond1:	# delete all commas so that times written as 2,454,000 -> 2454000 that regex recongizes as a number
				data = data.replace(",", "")
				numbers = [float(s) for s in re.findall(r'-?\d+\.?\d*', data)]
				#print('numbers: ', numbers)

				for i in range(len(numbers)):
					measurement = numbers[i]
					#if measurement > 2000000:
					#	print('Mid-time: ', measurement)
					if measurement > 2000000: #and (data.find("BJD") or data.find("HJD")) and (data.find("epoch") or data.find("mid-transit time") or data.find("transit time") or data.find("mid-point")):
						count +=1
						tex_files +=1
						IDs.append(arxiv_id)
						print('arxiv id containing number > 2000000: ', str(arxiv_id))
		
						try:
							paper = arxiv.query(id_list=[arxiv_id])[0]
							arxiv.download(paper, dirpath = direct + f'/arxiv/article_database/{planet}')
						except:
							paper = arxiv.query(id_list=['astro-ph/' + arxiv_id])[0]
							arxiv.download(paper, dirpath = direct + f'/arxiv/article_database/{planet}')

						break
	total_count +=1
	print(f'{arxiv_id}: # of tex files containing time > 2e6: ', tex_files)

            

print('count: ', count)
print('fraction from original number of articles: ', count/total_count)
np.savetxt(direct + f'/article_database/arxiv_ids_w_mid_times.txt', np.array(IDs), fmt='%s')


 