from __future__ import unicode_literals
# -*- coding: utf-8 -*-

__author__ = 'memery'

from bs4 import BeautifulSoup
import urllib2
from Bio import Entrez
import time
import re
import unicodecsv as csv


class Newsletter(object):

    def __init__(self, url):
        self.articles = []
        self.url = url
        self.soup = self.make_soup()
        self.publication_titles = self.parse_connexon()
        self.pubmed_list = self.look_up_titles()
        self.records = self.fetch_from_pubmed()
        self.parse_pubmed_soup()
        self.info = {'Lead Source': 'Connexon',
                     'Specific Lead Source': self.find_specific_lead_source(),
                     'Newsletter Archived Link': self.url.lstrip('http://www.'),
                     'Search Term': 'Connexon; {}'.format(self.find_specific_lead_source())}

    def make_soup(self):
        """Given a URL will return a BeatifulSoup of that URL

        Utilizes a header to avoid 503 Errors
        """
        header = {'User-Agent': "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/534.30 (KHTML, like Gecko) Ubuntu/11.04 "
                                "Chromium/12.0.742.112 Chrome/12.0.742.112 Safari/534.30"}
        request = urllib2.Request(self.url, headers=header)
        return BeautifulSoup(urllib2.urlopen(request))

    def parse_connexon(self):
        """Given a BeautifulSoup object, returns a list of publication names"""
        pubs = self.soup.find_all(_find_comment)
        publication_titles = []
        for pub in pubs:
            publication_titles.append(pub.text.lstrip('\n'))
        return publication_titles

    def look_up_titles(self):
        """Returns a list of PMIDs given a list of Connexon Titles"""
        pubmed_list = []
        for publication in self.publication_titles:
                pubmed_list.append(lookup_up_title(publication))

        return pubmed_list

    def fetch_from_pubmed(self):
        """Returns a list of Pubmed records based on a list of PMIDs"""
        Entrez.email = "matthew.emery@stemcell.com"
        handle = Entrez.efetch(db='pubmed', id=self.pubmed_list, retmode='xml')
        return BeautifulSoup(handle.read())

    def find_specific_lead_source(self):
        """Returns the name of the specific lead source."""
        title = self.soup.find('title').text
        title_split = title.split(' - ')
        vol = re.findall('\d+.\d+', title_split[0])
        return '{} {}'.format(title_split[1], vol[0])

    def parse_pubmed_soup(self):
        """Appends to self.list_of_articles Article objects"""
        for article in self.records('pubmedarticle'):
            self.articles.append((Article(article)))

    def write_csv(self, csv_file):
        """Adds line to a CSV contain all the information contained in self.articles"""
        field_names = ('First Name', 'Last Name', 'Email', 'Company', 'Department', 'Lead Source',
                       'Specific Lead Source', 'Newsletter Archived Link', 'Search Term', 'Publication Date',
                       'Publication Link', 'Article Title', 'Aff')
        csv_writer = csv.DictWriter(csv_file, fieldnames=field_names)
        csv_writer.writeheader()
        for article in self.articles:
            for author in article.authors:
                full_dict = dict(author.info.items() + article.info.items() + self.info.items())
                csv_writer.writerow(full_dict)

        csv_file.close()


class Article(object):

    def __init__(self, tag):
        self.tag = tag
        self.info = {}
        self.authors = []
        self.find_title()
        self.find_date()
        self.find_doi()
        self.find_authors()

    def find_title(self):
        self.info['Article Title'] = self.tag.articletitle.text.strip().strip('.')

    def find_date(self):
        year, month, day = None, None, None
        if self.tag.find('pubmedpubdate', {'pubstatus': 'aheadofprint'}):
            year = self.tag.find(pubstatus='aheadofprint').year.text.strip()
            month = self.tag.find(pubstatus='aheadofprint').month.text.strip()
            day = self.tag.find(pubstatus='aheadofprint').day.text.strip()
        elif self.tag.pubdate:
            try:
                year = self.tag.pubdate.year.text.strip()
                month = self.tag.pubdate.month.text.strip()
                day = self.tag.pubdate.day.text.strip()
            except AttributeError:
                year = self.tag.find(pubstatus='medline').year.text.strip()
                month = self.tag.find(pubstatus='medline').month.text.strip()
                day = self.tag.find(pubstatus='medline').day.text.strip()
        if month.isdigit():
            pubdate = time.strptime('{}{}{}'.format(month, day, year), '%m%d%Y')
        else:
            pubdate = time.strptime('{}{}{}'.format(month, day, year), '%b%d%Y')
        self.info['Publication Date'] = time.strftime('%m/%d/%Y', pubdate)

    def find_doi(self):
        """Return the DOI of an article as a string"""
        if self.tag.find(idtype='doi'):
            try:
                article_url = urllib2.urlopen('http://dx.doi.org/{}'.format(self.tag.find(idtype='doi').text.strip()))
                self.info['Publication Link'] = article_url.url
            except urllib2.HTTPError:
                self.info['Publication Link'] = 'DOI cannot be resolved: ' \
                                                'http://dx.doi.org/{}'.format(self.tag.find(idtype='doi').text.strip())
        else:
            self.info['Publication Link'] = 'DOI not found'

    def find_authors(self):
        prev_aff = ''

        for author in self.tag('author'):
            obj_author = Author(author)
            if obj_author.info['Aff'] == '':
                obj_author.set_institute(prev_aff)
            else:
                prev_aff = obj_author.info['Aff']
            self.authors.append(obj_author)


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
        self.info['Email'] = regex_search(self.info['Aff'], 'Email')

    def set_institute(self, aff):
        self.info['Aff'] = aff
        self.find_company()
        self.find_department()
        self.find_email()


def regex_search(institute, mode, lastname=''):
    """Attempt to parse the institute entry using regular expressions"""
    regex_dict = {'Department': r'[\w ]*Department[\w ]*|[\w ]*Laboratory[A-Z ]*|'
                                r'[\w ]*Cent[er|re][\w ]*|[\w ]*Service[A-Z ]*|[\w ]*Service[A-Z ]*'
                                r'|[\w ]*Institute[A-Z ]*',
                  'Company': r'[\w ]*Universit[y|aria][\w ]*|[\w \']*Institut[e]?[\w \']*|'
                             r'[\w ]*ETH[\w ]*|[\w \']*Academy[\w \']*|[\w \'&]*College[\w \']*',
                  'Email': r'\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,4}'
                  }
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


def lookup_up_title(publication_title):
    """Returns Entrez entry for a search of the publication title. If publication title does not return result,
    input PMID manually to continue to retrieve entry"""
    Entrez.email = "matthew.emery@stemcell.com"
    handle = Entrez.esearch(db='pubmed', term=publication_title, retmax=1)
    try:
        return Entrez.read(handle)['IdList'][0]
    except IndexError:
        # return manual_pmid(publication_title)
        # uncomment for actual running, I don't understand how mock works
        return '25430711'


def url_wrapper():
    """Prompts user for Connexon issue URL. Will raise AssertionError if non-standard URL added.

    Returns URL"""
    url = raw_input('Paste the URL of a Connexon website here: ')
    assert '/issue/' in url, 'URL does not point to Connexon issue.'
    return url

if __name__ == '__main__':
    url = url_wrapper()
    news = Newsletter(url)
    news.write_csv(open('leadentry.csv', 'wb'))

#TODO: Catch Press Release first titles
