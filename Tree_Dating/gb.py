"""
Usage:

import gb
entry = gb.GenbankEntry("JQ004916")
date = entry.source_dict["collection_date"]
"""

import sys
from time import sleep

if sys.version_info[0] == 3:
	from urllib.request import urlopen
else:
	from urllib import urlopen

from bs4 import BeautifulSoup

class GenbankEntry:
	
	def __init__(self, accno):
		self.accno = accno
		self.xml = pull_genbank_xml(accno)
		if self.xml is None:
			return	
		"""
		while self.xml is None:
			print ("Failed to pull xml for sequence %s: trying again." % accno)
			self.xml = pull_genbank_xml(accno)
		"""
		
		defn = self.xml.find("GBSeq_definition")

		try:
			self.name = self.xml.find("GBSeq_definition").string
		except AttributeError:
			return None
		except:
			raise
		
		self.pmids = []
		
		refs = self.xml.find_all("GBReference")
		for ref in refs:
		    pmid = ref.find('GBReference_pubmed')
		    if pmid is not None:
		        self.pmids.append(pmid.text)
		
		comment = self.xml.find("GBSeq_comment")
		if comment:
			self.comment = comment.string
		else:
			self.comment = ""

		self.source_dict = {}
		source = self.xml.find("GBFeature")
		for qual in source.find_all("GBQualifier"):
			key = qual.GBQualifier_name
			value = qual.GBQualifier_value
			if key is not None and value is not None:
				self.source_dict[key.string] = value.string

		self.sequence = self.xml.find("GBSeq_sequence").string

		self.genes = {}
		for feature in self.xml.find_all("GBFeature"):
			if feature.GBFeature_key.string == "gene":
				start = int(feature.GBFeature_intervals.GBInterval.GBInterval_from.string)-1
				end = int(feature.GBFeature_intervals.GBInterval.GBInterval_to.string)
				for qualifier in feature.find_all("GBQualifier"):
					if qualifier.GBQualifier_name.string == "gene":
						gene = qualifier.GBQualifier_value.string
						self.genes[gene] = self.sequence[start:end]
			
	def add_products(self, db):
		for gene in self.genes:
			db.add_product(self.accno, gene, self.genes[gene])

def pull_pubmed_xml(pmid):
	return pull_xml("pubmed", "null", "xml", pmid)

def pull_genbank_xml(accno):
	# First pull the genbank ID for the sequnence (which is not the same as the
	# accession number, for some reason).
	gb_url = "http://www.ncbi.nlm.nih.gov/nuccore/{0}".format(accno)
	gb_html = None
	while gb_html is None:
		try:
			gb_html = BeautifulSoup(urlopen(gb_url), 'html.parser')
		except IOError:
			print ('HTML protocol error, retrying after 10 seconds')
			sleep(10)
			continue
		except:
			raise
	
	gbid_tag = gb_html.find(id="viewercontent1")
	if gbid_tag is not None:
		gbid = gbid_tag.get("val")
		return pull_xml("nuccore", "gb", "xml", gbid)

def pull_xml(db, rettype, retmode, id):
	base_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
	url_template = "{0}?db={1}&rettype={2}&retmode={3}&id={4}"
	
	url = url_template.format(base_url, db, rettype, retmode, id)
	return BeautifulSoup(urlopen(url), "xml")
