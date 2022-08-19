# pip3 install requests beautifulsoup4

import re
import time
import json
import requests
from bs4 import BeautifulSoup

N = 20  # number of problems
PARAMS = {
	'DNA': ['s'],
	'RNA': ['t'],
	'REVC': ['s'],
	'FIB': ['n', 'k'],
	'GC': ['fasta'],
	'HAMM': ['s', 't'],
	'IPRB': ['k', 'm', 'n'],
	'PROT': ['s'],
	'SUBS': ['s', 't'],
	'CONS': ['fasta'],
	'FIBD': ['n', 'm'],
	'GRPH': ['fasta'],
	'IEV': ['a', 'b', 'c', 'd', 'e', 'f'],
	'LCSM': ['fasta'],
	'LIA': ['k', 'N'],
	'MPRT': ['ids'],
	'MRNA': ['protein'],
	'ORF': ['s'],
	'PERM': ['n'],
	'PRTM': ['P'],
}


# ------------------------------------------------------------------------------
# Get list of problems (Bioinformatics Stronghold)
# ------------------------------------------------------------------------------

r = requests.get("http://rosalind.info/problems/list-view/")
soup = BeautifulSoup(r.content, "html.parser")
soup_problems = soup.find_all("tr")[1:]
problems = []
for problem in soup_problems:
	problems.append({
		"id": problem.find_all("td")[0].get_text(),
		"title": problem.find_all("td")[1].find("a").get_text().strip(),
		# "url": "http://rosalind.info" + problem.find_all("td")[1].find("a").get("href")
	})

# ------------------------------------------------------------------------------
# Get each problem's content
# ------------------------------------------------------------------------------

for problem in problems[0:N]:
	r = requests.get(f"http://rosalind.info/problems/{problem['id'].lower()}")
	time.sleep(0.2)  # don't overwhelm the server
	soup = BeautifulSoup(r.content, "html.parser")

	# Parse given/return
	problem['given'] = mathjax_to_md(soup.select(".given-return")[0].parent.get_text().strip()).replace('Given: ', '').replace('\n', ' ')
	problem['return'] = mathjax_to_md(soup.select(".given-return")[1].parent.get_text().strip()).replace('Return: ', '').replace('\n', ' ')
	# Parse sample input/output
	problem['sample_data'] = soup.select(".codehilite pre")[0].get_text().strip()
	problem['sample_output'] = soup.select(".codehilite pre")[1].get_text().strip()

	# Track function parameters
	problem['params'] = PARAMS[problem['id']]


# ------------------------------------------------------------------------------
# Output result
# ------------------------------------------------------------------------------

# Manual fixes
problems[4]['return'] = problems[4]['return'].replace('; please see the note on absolute error below', '')  # note isn't there :)
problems[5]['given'] = problems[5]['return'] = 'The Hamming distance _dH(s, t)_.'  # remove extra italic
problems[11]['return'] = problems[11]['return'].replace('_O_3_', '_O3_')
problems[15]['sample_data'] = 'A2Z669\nB5ZC00\nP07204_TRBM_HUMAN\nP20840_SAG1_YEAST'
problems[15]['sample_output'] = 'B5ZC00\n85 118 142 306 395\nP07204_TRBM_HUMAN\n47 115 116 382 409\nP20840_SAG1_YEAST\n79 109 135 248 306 348 364 402 485 501 614'

print(json.dumps(problems[0:N]).replace('{', '\n{'))


# ------------------------------------------------------------------------------
# Utilities
# ------------------------------------------------------------------------------

def mathjax_to_md(string):
	for group in re.findall('\$.+?\$', string):
		string = string.replace(group, "_" + group[1:-1] + "_")
	string = re.sub('\{\\\\textrm\{(.+?)\}\}', r'\1', string)
	string = re.sub('\{\\\\mathrm\{(.+?)\}\}', r'\1', string)
	string = re.sub('\\\\mathrm\{(.+?)\}', r'\1', string)
	string = re.sub('\\\\mathscr\{(.+?)\}', r'\1', string)
	string = re.sub('\\\\ldots', '...', string)
	string = re.sub('\\\\{', '{', string)
	string = re.sub('\\\\}', '}', string)
	string = string.replace('\\leq', '<=')
	return string
