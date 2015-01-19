#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import unicode_literals

__author__ = 'memery'

from bs4 import BeautifulSoup
import urllib2
from Bio import Entrez
import time
import re
import unicodecsv as csv
import socket
import Tkinter
import tkFileDialog


class Batch(object):
    """A Batch object contains a collection of PMIDs that are parsed into Articles"""

    def __init__(self):
        self.articles = []
        self.pubmed_xml = None

    def add_article(self, article):
        self.articles.append(article)

    def create_pubmed_xml(self):
        """Returns a BeautifulSoup object from a list of Pubmed IDs.

        Creates one BeautifulSoup object from a list of Pubmed IDs. The Entrez email is currently defaulted to
        matthew.emery@stemcell.com. {This may change in the future}. The Beautiful Soup is retrieved in XML format.
        I believe this may affect whether Greek letters are properly encoded.
        """
        Entrez.email = "matthew.emery@stemcell.com"
        queries = [article['PMID'] for article in self.articles if 'PMID' in article]
        handle = Entrez.efetch(db='pubmed', id=queries, retmode='xml')
        self.pubmed_xml = BeautifulSoup(handle.read())

    def parse_pubmed_soup(self):
        """Appends to self.list_of_articles Article objects"""
        if self.pubmed_xml:
            for article in self.articles:
                if article.in_info('PMID'):
                    child = self.pubmed_xml.find('pmid', text=re.compile(article.get_info('PMID')))
                    article.update_info_dict('Tag', child.find_parent('pubmedarticle'))
        else:
            raise AssertionError('Can\'t parse a PubMed soup that isn\'t there. Try self.create_pubmed_xml')


class ZoteroEntry(Batch):
    def __init__(self, zotero_csv, info=None):
        Batch.__init__(self)
        self.zotero_csv = zotero_csv
        if info:
            self.info = info
        else:
            self.info = {'Lead Source': 'Web Search (Google, FASEB, PubMed, CRISP)'}

    def read_csv(self):
        with self.zotero_csv as zotero:
            reader = csv.DictReader(zotero)
            for row in reader:
                if row['Item Type'] == 'journalArticle':
                    article = Article(info={'Article Title': row['Title'],
                                            'Publication Link': row['Url']})
                    article.update_info_dict('PMID', article.lookup_up_pmid())
                    if not article.in_info('PMID'):
                        article.update_info_dict('Authors', row['Author'].split('; '))
                    self.add_article(article)

    def write_csv(self, csv_file):
        """Complete"""
        field_names = ('Publication Link', 'Publication Date', 'Article Title', 'Abstract', 'Search Term',
                       'Product Use/Assay Type', 'Area of Interest', 'Product Line', 'First Name', 'Last Name',
                       'Company', 'Department', 'Email', 'Lead Source', 'Specific Lead Source', 'Product Sector',
                       'PMID')

        csv_writer = csv.DictWriter(csv_file, extrasaction='ignore', fieldnames=field_names)
        csv_writer.writeheader()
        for article in self.articles:
            for author in article.authors:
                full_dict = dict(author.info.items() + article.info.items() + self.info.items())
                csv_writer.writerow(full_dict)

        csv_file.close()


class Newsletter(Batch):

    def __init__(self, url):
        Batch.__init__(self)
        self.url = url
        self.soup = self.make_soup()
        self.publication_titles = self.parse_connexon()
        for pub in self.publication_titles:
            self.lookup_up_title(pub)
        self.pubmed_xml = self.create_pubmed_xml()
        print self.pmids
        self.parse_pubmed_soup()
        self.info = {'Lead Source': 'Connexon',
                     'Specific Lead Source': self.find_specific_lead_source(),
                     'Newsletter Archived Link': self.url.lstrip('http://www.'),
                     'Search Term': 'Connexon; {}'.format(self.find_specific_lead_source())}

    def make_soup(self):
        """Given a URL will return a BeautifulSoup of that URL

        Utilizes a header to avoid 503 Errors
        """
        header = {'User-Agent': "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/534.30 (KHTML, like Gecko) Ubuntu/11.04 "
                                "Chromium/12.0.742.112 Chrome/12.0.742.112 Safari/534.30"}
        request = urllib2.Request(self.url, headers=header)
        return BeautifulSoup(urllib2.urlopen(request, timeout=20))

    def parse_connexon(self):
        """Given a BeautifulSoup object, returns a list of publication names"""
        pubs = self.soup.find_all(_find_comment)
        publication_titles = []
        for pub in pubs:
            publication_titles.append(pub.text.lstrip('\n'))
        return publication_titles

    def find_specific_lead_source(self):
        """Returns the name of the specific lead source."""
        title = self.soup.find('title').text
        title_split = title.split(' - ')
        vol = re.findall('\d+.\d+', title_split[0])
        return '{} {}'.format(title_split[1], vol[0])

    def write_csv(self, csv_file):
        """Adds line to a CSV contain all the information contained in self.articles"""
        field_names = ('First Name', 'Last Name', 'Email', 'Company', 'Department', 'Lead Source',
                       'Specific Lead Source', 'Newsletter Archived Link', 'Search Term', 'Publication Date',
                       'Publication Link', 'Article Title', 'Aff', 'PMID')
        csv_writer = csv.DictWriter(csv_file, extrasaction='ignore', fieldnames=field_names)
        csv_writer.writeheader()
        for article in self.articles:
            for author in article.authors:
                full_dict = dict(author.get_info_items() + article.get_info_items() + self.info.items())
                csv_writer.writerow(full_dict)

        csv_file.close()


class Article(object):

    def __init__(self, info=None):

        # self.tag = tag
        if info:
            self.info = info
        else:
            self.info = {}
        # self.find_title()
        # self.find_date()
        # self.find_doi()
        # self.find_authors()
        # self.find_abstract()

    def in_info(self, query):
        return query in self.info

    def get_info(self, query):
        if query in self.info:
            return self.info[query]
        else:
            return None

    def lookup_up_pmid(self, translated=False):
        """Returns Entrez entry for a search of the publication title. If publication title does not return result,
        input PMID manually to continue to retrieve entry"""
        if self.info['Article Title']:
            Entrez.email = "matthew.emery@stemcell.com"
            handle = Entrez.esearch(db='pubmed', term=self.info['Article Title'], retmax=10, sort='relevance')
            try:
                return Entrez.read(handle)['IdList'][0]
            except IndexError:
                if not translated:
                    # Britishness could be the issue
                    brit_dict = {'Leukemia': 'Leukaemia',
                                 'Tumor': 'Tumour',
                                 'Signaling': 'Signalling',
                                 'α': 'alpha'}
                    for brit in brit_dict.items():
                        publication_title = self.info['Article Title'].replace(brit[0], brit[1])
                    self.update_info_dict('Article Title', publication_title)
                    self.lookup_up_pmid(self, publication_title, translated=True)
                else:
                    print 'Could not find PMID: {}'.format(self.info['Article Title'])
        else:
            raise AssertionError('Impossible to search PubMed without a title!')

    def update_info_dict(self, key, value):
        self.info[key] = value

    def get_info_items(self):
        return self.info.items()

    def process_tag_info(self, tag_list):
        """Looks for tag info for a list of tags"""

    def find_title(self):
        self.info['Article Title'] = self.info['Tag'].articletitle.text.strip().strip('.')

    def find_date(self): #This sucks. Rewrite it
        year, month, day = None, None, None
        if self.info['Tag'].find('pubmedpubdate', {'pubstatus': 'aheadofprint'}):
            year = self.info['Tag'].find(pubstatus='aheadofprint').year.text.strip()
            month = self.info['Tag'].find(pubstatus='aheadofprint').month.text.strip()
            day = self.info['Tag'].find(pubstatus='aheadofprint').day.text.strip()
        elif self.info['Tag'].pubdate:
            try:
                year = self.info['Tag'].pubdate.year.text.strip()
                month = self.info['Tag'].pubdate.month.text.strip()
                day = self.info['Tag'].pubdate.day.text.strip()
            except AttributeError:
                year = self.info['Tag'].find(pubstatus='medline').year.text.strip()
                month = self.info['Tag'].find(pubstatus='medline').month.text.strip()
                day = self.info['Tag'].find(pubstatus='medline').day.text.strip()
        if month.isdigit():
            pubdate = time.strptime('{}{}{}'.format(month, day, year), '%m%d%Y')
        else:
            pubdate = time.strptime('{}{}{}'.format(month, day, year), '%b%d%Y')
        self.info['Publication Date'] = time.strftime('%m/%d/%Y', pubdate)

    def find_doi(self):
        """Return the DOI of an article as a string"""
        if self.info['Tag'].find(idtype='doi'):
            try:
                doi = self.info['Tag'].find(idtype='doi').text.strip()
                header = {'User-Agent': "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/534.30 (KHTML, like Gecko) "
                                        "Ubuntu/11.04 Chromium/12.0.742.112 Chrome/12.0.742.112 Safari/534.30"}
                request = urllib2.Request('http://dx.doi.org/{}'.format(doi), headers=header)
                article_url = urllib2.urlopen(request, timeout=20)
                self.info['Publication Link'] = article_url.geturl()
            except (urllib2.HTTPError, socket.timeout) as e:
                self.info['Publication Link'] = '{}: http://dx.doi.org/{}'.format(e, self.info['Tag'].find(idtype='doi').text)
        else:
            self.info['Publication Link'] = 'DOI not found'

    def find_authors(self):
        prev_aff = ''

        for author in self.info['Tag']('author'):
            obj_author = Author(author)
            if obj_author.info['Aff'] == '':
                obj_author.set_institute(prev_aff)
            else:
                prev_aff = obj_author.info['Aff']
            self.authors.append(obj_author)
            print obj_author

    def find_abstract(self):
        self.info['Abstract'] = self.tag.abstracttext.text.strip()

    def __str__(self):
        return self.info['Article Title'].encode('UTF-8')

    def __repr__(self):
        return self.info['Article Title'].encode('UTF-8')


class Author(object):

    def __init__(self, tag):
        self.tag = tag
        self.info = {}
        self.find_first_name()
        self.find_last_name()
        self.find_institute()
        self.find_department()
        self.find_company()
        self.find_email()

    def find_last_name(self):
        self.info['Last Name'] = unicode(self.tag.lastname.text.strip())

    def find_first_name(self):
        forename = self.tag.forename.text.strip()
        if forename.split(' '):
            self.info['First Name'] = unicode(forename.split(' ')[0])
        else:
            self.info['First Name'] = unicode(forename)

    def find_institute(self):
        try:
            self.info['Aff'] = self.tag.affiliation.text.strip()
        except AttributeError:
            self.info['Aff'] = ''

    def find_department(self):
        self.info['Department'] = regex_search(self.info['Aff'], 'Department')

    def find_company(self):
        self.info['Company'] = regex_search(self.info['Aff'], 'Company')

    def find_email(self):
        self.info['Email'] = regex_search(self.info['Aff'], 'Email', lastname=self.info['Last Name'])

    def set_institute(self, aff):
        self.info['Aff'] = aff
        self.find_company()
        self.find_department()
        self.find_email()

    def get_info_items(self):
        return self.info.items()

    def __str__(self):
        return '{} {}'.format(self.info['First Name'], self.info['Last Name']).encode('UTF-8')


def regex_search(institute, mode, lastname=''):
    """Attempt to parse the institute entry using regular expressions"""
    regex_dict = {'Department': r'[\w ]*Department[\w ]*|[\w ]*Laboratory[A-Z ]*|'
                                r'[\w ]*Cent[er|re][\w ]*|[\w ]*Service[A-Z ]*|[\w ]*Service[A-Z ]*'
                                r'|[\w ]*Institute[A-Z ]*',
                  'Company': r'[\w\- ]*Universit[y|aria][\w ]*|[\w \']*Institut[e]?[\w \']*|'
                             r'[\w ]*ETH[\w ]*|[\w \']*Academy[\w \']*|[\w \'&]*College[\w \']*',
                  'Email': r'\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,4}'}
    query = re.findall(regex_dict.get(mode), institute, flags=re.I | re.U)
    if query:
        if mode == 'Department':
            return query[0].strip()
        elif mode == 'Company':
            return query[-1].strip()
        else:
            if len(query) == 1:
                return query[0]
            else:
                for email in query:
                    if lastname.lower() in email.lower():
                        return email
                return str(query)
    else:
        return ''


def _find_comment(tag):
    """Don't use this function directly. Is used in find_comment to pull all lines of text with #PUBLICATION TITLE
    comment"""
    try:
        return tag.has_attr('face') and tag.has_attr('size') and '#PUBLICATIONS TITLE' in tag.contents[1]
    except IndexError:
        return False

# def _find_pmid_in_xml(tag):
#     try:



def url_wrapper():
    """Prompts user for Connexon issue URL. Will raise AssertionError if non-standard URL added.

    Returns URL"""
    url = raw_input('Paste the URL of a Connexon website here: ')
    assert '/issue/' in url, 'URL does not point to Connexon issue.'
    return url


def make_zotero_entry():
    info_dict = {'Search Term': raw_input('Input Search Term: '),
                 'Product Use/Assay Type': raw_input('Input Product Use/Assay Type: '),
                 'Product Line': raw_input('Input Product Line: '),
                 'Area of Interest': raw_input('Input Area of Interest: '),
                 'Product Sector': raw_input('Input Product Sector: '),
                 'Lead Source': 'Web Search (Google, FASEB, PubMed, CRISP)'}
    root = Tkinter.Tk()
    root.withdraw()
    print 'Please Select Zotero CSV.'
    source = tkFileDialog.askopenfile(parent=root,
                                      title='Select Zotero CSV')
    return ZoteroEntry(source, info_dict)


if __name__ == '__main__':
    PROMPT = raw_input('Press "Z" for Zotero. Press "C" for Connexon. ')
    if PROMPT.upper() == 'Z':
        batch = make_zotero_entry()
        batch.read_csv()
    elif PROMPT.upper() == 'C':
        chosen_url = url_wrapper()
        batch = Newsletter(chosen_url)
    else:
        raise AssertionError('Invalid choice, try again')
    batch.pubmed_xml = batch.create_pubmed_xml()
    batch.parse_pubmed_soup()
    batch.write_csv(open('BatchOutput.csv', 'wb'))


    # TODO: Catch Press Release first titles
    # TODO: Search for Salesforce IDs?
    # TODO: Custom Entrez.email entries
    # TODO: A reminder of where you are: self.pubmed_xml.find('pmid', text=re.compile('25433608')).parent.find('issn').text