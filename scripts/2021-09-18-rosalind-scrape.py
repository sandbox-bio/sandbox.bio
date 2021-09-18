# pip3 install requests beautifulsoup4

import re
import time
import requests
from bs4 import BeautifulSoup

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
		"url": "http://rosalind.info" + problem.find_all("td")[1].find("a").get("href")
	})

# ------------------------------------------------------------------------------
# Get each problem's content
# ------------------------------------------------------------------------------

for problem in problems:
	r = requests.get(problem["url"])
	time.sleep(0.5)
	soup = BeautifulSoup(r.content, "html.parser")

	problem['given'] = mathjax_to_md(soup.select(".given-return")[0].parent.get_text()).replace('Given: ', '')
	problem['return'] = mathjax_to_md(soup.select(".given-return")[1].parent.get_text()).replace('Return: ', '')

	problem['sample_data'] = soup.select(".codehilite pre")[0].get_text()
	problem['sample_output'] = soup.select(".codehilite pre")[1].get_text()


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
