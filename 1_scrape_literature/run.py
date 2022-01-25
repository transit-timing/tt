import shutil
import argparse
import json
import os
import re
import sys
import tarfile
import tarfile
import lxml.html
import requests
import numpy as np


com_pt = re.compile(r'(?<!\\)%+(.+)')
multi_com_pt = re.compile(r'\\begin{comment}(.+?)\\end{comment}')

arxiv_id_pt = re.compile(r'(?<!\d)(\d{4}\.\d{5})(?!\d)')
url_base = 'https://arxiv.org/e-print/'
url_bib_base = 'http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:'

direct = os.path.dirname(os.getcwd())

def get_all_arxiv_ids(text):
    ids = []
    for arxiv_id in arxiv_id_pt.findall(text):
        ids.append(arxiv_id)
    return ids

def download(url, dir_path='./'):
    idx = os.path.split(url)[-1]
    file_name = idx + '.tar.gz'
    file_path = os.path.join(dir_path, file_name)
    if os.path.exists(file_path):
        return file_path

    r = requests.get(url)
    sys.stderr.write('\tdownload {}'.format(url) + '\n')
    if r.status_code == 200:
        with open(file_path, 'wb') as f:
            f.write(r.content)
        return file_path
    else:
        return 0

def read_papers(arxiv_ids, dir_path='./'):
    results = {}
    for arxiv_id in arxiv_ids:
        sys.stderr.write('[{}]'.format(arxiv_id) + '\n')
        result = read_paper(arxiv_id, dir_path)
        if result:
            if 'title' in result:
                sys.stderr.write('\t({})'.format(result['title']) + '\n')
                sys.stderr.write('\t {}'.format(' / '.join(result['authors'])) + '\n')
            results[arxiv_id] = result
    return results

def read_paper(arxiv_id, dir_path='./'):
    #print('url base: ', url_base)
    url = url_base + arxiv_id
    targz_path = download(url, dir_path)

    #url = url_bib_base + arxiv_id
    #print('url base: ', url_bib_base)
    #bib_path = download(url, dir_path)



def extract_bibinfo(f):
    info = {'title': '', 'authors': []}
    dom = lxml.html.fromstring(f.read())
    for c in dom.xpath('//meta'):
        name = c.attrib.get('name', '')
        if name == 'dc.title':
            info['title'] = c.attrib['content']
        elif name == 'dc.creator':
            info['authors'].append(c.attrib['content'])
    return info

def extract_comment(f):
    results = []
    for line_idx, line in enumerate(f.readlines()):
        for comment in com_pt.findall(line.decode('utf-8')):
            results.append(comment)

    for comment in multi_com_pt.findall(f.read().decode('utf-8')):
        results.append(comment)

    return results

def main():
    parser = argparse.ArgumentParser(description='Arxiv')
    parser.add_argument('--text', '-t', help='text which contains arxiv ids', default = direct + '/arxiv/article_database/arxiv_ids.txt')
    parser.add_argument('--id', '-i', nargs='+', default=[])
    parser.add_argument('--save-dir', '-s', default='./')
    parser.add_argument('--output', '-o', default='./comments.json')
    parser.add_argument('--planet', '-p')

    

    args = parser.parse_args()
    sys.stderr.write(json.dumps(args.__dict__, indent=2) + '\n')

    planet_name = args.planet

    save_dir = direct + f'/arxiv/arxiv_dir/{planet_name}'

    if not os.path.isdir(save_dir): 
        os.mkdir(save_dir)


    if not os.path.isdir(direct + f'/arxiv/untar_dir/{planet_name}/'): 
        os.mkdir(direct + f'/arxiv/untar_dir/{planet_name}/')


    ids = args.id
    if args.text:
        sys.stderr.write('load text: {}'.format(args.text) + '\n')
        #ids.extend(get_all_arxiv_ids(open(args.text).read()))
    #ids = np.loadtxt(args.text)
    ids = np.genfromtxt(args.text,dtype='str')
    #ids = np.genfromtxt(direct + '/arxiv/article_database/arxiv_ids.txt',dtype='str')
    #print('IDs: ', ids)

    #sys.stderr.write('TARGET:\n' + '\n'.join('{} {}'.format(i, idx) for i, idx in enumerate(ids)) + '\n\n')
    all_results = read_papers(ids, save_dir)

    #print(json.dumps(all_results, indent=2))
    #json.dump(all_results, open(args.output, 'w'), indent=2)
    untarable = []

    for filename in os.listdir(direct + f'/arxiv/arxiv_dir/{planet_name}'):
        
        if filename.endswith(".tar.gz"):
            #print('filename: ', filename)

            try:

                tar = tarfile.open(direct + f'/arxiv/arxiv_dir/{planet_name}/' + filename, "r:gz")

                name = os.path.splitext(filename)[0]
                name = os.path.splitext(name)[0]

                if not os.path.isdir(direct + f'/arxiv/untar_dir/{planet_name}/'): 
                    os.mkdir(direct + f'/arxiv/untar_dir/{planet_name}/')

                    if not os.path.isdir(direct + f'/arxiv/untar_dir/{planet_name}/{name}/'): 
                        os.mkdir(direct + f'/arxiv/untar_dir/{planet_name}/{name}/')


                tar.extractall(path=direct + f'/arxiv/untar_dir/{planet_name}/{name}/')
                tar.close()
            except:
                print('Could not unzip: ', name)



def extract_all(archives, extract_path):
    for filename in archives:
        shutil.unpack_archive(filename, extract_path)

if __name__ == '__main__':
    main()

