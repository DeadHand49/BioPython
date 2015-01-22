#!/usr/bin/python
# -*- coding: utf-8 -*-
"""A script that uses PubMed to fill in lead entry on a csv."""

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

    def __init__(self, info=None, field_names=None):
        self.articles = []
        self.pubmed_xml = None
        if info:
            self.info = info
        else:
            self.info = {'Object Type': 'Batch Test'}
        if field_names:
            self.field_names = field_names
        else:
            self.field_names = ('Article Title', 'PMID', 'Last Name', 'First Name', 'Email', 'Company', 'Department')

    def add_article(self, article):
        assert isinstance(article, Article), 'Only Articles go in self.articles'
        print article
        self.articles.append(article)

    def create_pubmed_xml(self):
        """Returns a BeautifulSoup object from a list of Pubmed IDs.

        Creates one BeautifulSoup object from a list of Pubmed IDs. The Entrez email is currently defaulted to
        matthew.emery@stemcell.com. {This may change in the future}. The Beautiful Soup is retrieved in XML format.
        I believe this may affect whether Greek letters are properly encoded.
        """
        Entrez.email = "matthew.emery@stemcell.com"
        queries = [article.get_info('PMID') for article in self.articles if article.in_info('PMID')]
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

    def write_csv(self, csv_file):
        """Complete"""
        print 'Writing CSV.'
        csv_writer = csv.DictWriter(csv_file, extrasaction='ignore', fieldnames=self.field_names)
        csv_writer.writeheader()
        for article in self.articles:
            for author in article.authors:
                full_dict = dict(author.info.items() + article.info.items() + self.info.items())
                csv_writer.writerow(full_dict)

        csv_file.close()


class ZoteroEntry(Batch):
    def __init__(self, zotero_csv, info=None, field_names=None):
        Batch.__init__(self)
        self.zotero_csv = zotero_csv
        if info:
            self.info = info
        else:
            self.info = {'Lead Source': 'Web Search (Google, FASEB, PubMed, CRISP)'}
        if field_names:
            self.field_names = field_names
        else:
            self.field_names = ('Publication Link', 'Publication Date', 'Article Title', 'Abstract', 'Search Term',
                                'Product Use/Assay Type', 'Area of Interest', 'Product Line', 'First Name', 'Last Name',
                                'Company', 'Department', 'Email', 'Lead Source', 'Specific Lead Source',
                                'Product Sector', 'PMID')

    def read_csv(self):
        with self.zotero_csv as zotero:
            reader = csv.DictReader(zotero)
            for row in reader:
                if row['Item Type'] == 'journalArticle':
                    article = Article(info={'Article Title': row['Title'],
                                            'Publication Link': row['Url']})
                    article.update_info_dict('PMID', article.lookup_up_pmid(article.info['Article Title']))
                    if not article.in_info('PMID'):
                        article.update_info_dict('Authors', row['Author'].split('; '))
                    self.add_article(article)

    def construct_articles(self):
        for article in self.articles:
            if article.in_info('Tag'):
                article.update_info_dict('Publication Date', article.find_date())
                article.update_info_dict('Abstract', article.find_abstract())
                article.update_info_dict('Authors', article.find_authors())


class Newsletter(Batch):
    def __init__(self, url, info=None, field_names=None):
        Batch.__init__(self)
        self.url = url
        self.soup = self.make_soup()
        self.parse_connexon()
        if info:
            self.info = info
        else:
            self.info = {'Lead Source': 'Connexon',
                         'Specific Lead Source': self.find_specific_lead_source(),
                         'Newsletter Archived Link': self.url.lstrip('http://www.'),
                         'Search Term': 'Connexon; {}'.format(self.find_specific_lead_source())}
        if field_names:
            self.field_names = field_names
        else:
            self.field_names = ('First Name', 'Last Name', 'Email', 'Company', 'Department', 'Lead Source',
                                'Specific Lead Source', 'Newsletter Archived Link', 'Search Term', 'Publication Date',
                                'Publication Link', 'Article Title', 'Aff', 'PMID')

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
        pubs = self.soup.find_all(_find_comment)  # pubs[0].find_previous('a')
        for pub in pubs:
            article = Article(info={'Article Title': pub.text.lstrip('\n')})
            article.update_info_dict('PMID', article.lookup_up_pmid(pub.text.lstrip('\n')))
            article.update_info_dict('Publication Link', pub.find_previous('a').get('href'))
            self.add_article(article)

    def find_specific_lead_source(self):
        """Returns the name of the specific lead source."""
        title = self.soup.find('title').text
        return title.split(' - ')[1]

    def construct_articles(self):
        for article in self.articles:
            if article.in_info('Tag'):
                article.update_info_dict('Publication Date', article.find_date())
                article.update_info_dict('Authors', article.find_authors())


class Article(object):
    def __init__(self, info=None):
        if info:
            self.info = info
        else:
            self.info = {}
        self.authors = []

    def in_info(self, query):
        try:
            if self.info[query]:
                return query in self.info
            else:
                return False
        except KeyError:
            return False

    def get_info(self, query):
        if query in self.info:
            return self.info[query]
        else:
            return None

    def lookup_up_pmid(self, pub_title, translated=False):
        """Returns Entrez entry for a search of the publication title. If publication title does not return result,
        input PMID manually to continue to retrieve entry"""
        Entrez.email = "matthew.emery@stemcell.com"
        handle = Entrez.esearch(db='pubmed', term=pub_title, retmax=10, sort='relevance')
        try:
            print 'Found PMID: {}'.format(pub_title)
            return Entrez.read(handle)['IdList'][0]
        except IndexError:
            if not translated:
                pub_title = self.translate_british(pub_title)
                return self.lookup_up_pmid(pub_title, translated=True)
            else:
                print 'Could not find PMID: {}'.format(pub_title)

    @staticmethod
    def translate_british(publication_title):
        brit_dict = {'Leukemia': 'Leukaemia',
                     'Tumor': 'Tumour',
                     'Signaling': 'Signalling',
                     'α': 'alpha'}
        for brit in brit_dict.items():
            publication_title = publication_title.replace(brit[0], brit[1])
        return publication_title

    def update_info_dict(self, key, value):
        self.info[key] = value

    def get_info_items(self):
        return self.info.items()

    def find_title(self):
        return self.info['Tag'].articletitle.text.strip().strip('.')

    def find_date(self):
        pass
        # This is a mess figure it out later
        # year, month, day = None, None, None
        # potential_tags = [self.info['Tag'].find('pubmedpubdate', {'pubstatus': 'aheadofprint'}),
        # self.info['Tag'].pubdate,
        #                   self.info['Tag'].find(pubstatus='medline')]
        # potential_tags = [pot_tag for pot_tag in potential_tags if pot_tag]
        # while not day and not month and not year:
        #     for tag in potential_tags:
        #         if tag:
        #             year, month, day = self.return_date_from_tag(tag)
        #             if year and month and day:
        #                 break
        # return self.output_date(day, month, year)

    @staticmethod
    def output_date(day, month, year):
        if month.isdigit():
            pubdate = time.strptime('{}{}{}'.format(month, day, year), '%m%d%Y')
        else:
            pubdate = time.strptime('{}{}{}'.format(month, day, year), '%b%d%Y')
        return time.strftime('%m/%d/%Y', pubdate)

    @staticmethod
    def return_date_from_tag(find_tag):
        try:
            year = find_tag.year.strip()
            month = find_tag.month.strip()
            day = find_tag.day.strip()
            return year, month, day
        except (TypeError, AttributeError):  # consider a print statement here
            return None, None, None

    def find_doi(self):
        """Return the DOI of an article as a string"""
        if self.info['Tag'].find(idtype='doi'):
            doi = self.info['Tag'].find(idtype='doi').text.strip()
            try:
                header = {'User-Agent': "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/534.30 (KHTML, like Gecko) "
                                        "Ubuntu/11.04 Chromium/12.0.742.112 Chrome/12.0.742.112 Safari/534.30"}
                request = urllib2.Request('http://dx.doi.org/{}'.format(doi), headers=header)
                article_url = urllib2.urlopen(request, timeout=20)
                return article_url.geturl()
            except (urllib2.HTTPError, socket.timeout) as e:
                return '{}: http://dx.doi.org/{}'.format(e, doi)
        else:
            return 'DOI not found'

    def find_authors(self):
        prev_aff = ''

        for author_tag in self.info['Tag']('author'):
            author = Author(author_tag)
            author.update_info_dict('First Name', author.find_first_name())
            author.update_info_dict('Last Name', author.find_last_name())
            try:
                affiliation = author.find_affiliation()
                author.update_info_dict('Aff', affiliation)
                prev_aff = affiliation
            except AttributeError:
                author.update_info_dict('Aff', prev_aff)
            author.update_info_dict('Company', author.find_company())
            author.update_info_dict('Department', author.find_department())
            author.update_info_dict('Email', author.find_email())

            self.authors.append(author)
            print author

    def find_author(self):  # This needs to be filled in
        pass

    def find_abstract(self):
        return self.info['Tag'].abstracttext.text.strip()

    def __str__(self):
        return self.info['Article Title'].encode('UTF-8')

    def __repr__(self):
        return self.info['Article Title'].encode('UTF-8')


class Author(object):
    def __init__(self, tag):
        self.info = {'Tag': tag}

    def update_info_dict(self, key, value):
        self.info[key] = value

    def find_last_name(self):
        return unicode(self.info['Tag'].lastname.text.strip())

    def find_first_name(self):
        forename = self.info['Tag'].forename.text.strip()
        if forename.split(' '):
            return unicode(forename.split(' ')[0])
        else:
            return unicode(forename)

    def find_affiliation(self):
        return self.info['Tag'].affiliation.text.strip()

    def find_department(self):
        return regex_search(self.info['Aff'], 'Department')

    def find_company(self):
        return regex_search(self.info['Aff'], 'Company')

    def find_email(self):
        return regex_search(self.info['Aff'], 'Email', lastname=self.info['Last Name'])

    def in_info(self, query):
        return query in self.info

    def get_info_items(self):
        return self.info.items()

    def __str__(self):
        return '{} {}'.format(self.info['First Name'], self.info['Last Name']).encode('UTF-8')


def regex_search(institute, mode, lastname=''):  # possibly rafactor this
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
    batch.create_pubmed_xml()
    batch.parse_pubmed_soup()
    batch.construct_articles()
    batch.write_csv(open('BatchOutput.csv', 'wb'))


    # TODO: Catch Press Release first titles
    # TODO: Search for Salesforce IDs?
    # TODO: Custom Entrez.email entries