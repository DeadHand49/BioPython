#!/usr/bin/python
# -*- coding: utf-8 -*-
"""A script that uses PubMed to fill in lead entry on a csv.

This module allows the user to select either a previous CSV generated from Zotero or a Connexon URL as the basis for
the new leadentry CSV. The leadentry CSV will appear in the the same directory as the module with the title
'BatchOutput.csv'.

This module contains a simple structure. A Batch object (whether ZoteroEntry or Newsletter) acts as a factory,
building Article objects from previously known information and PubMed lookups. In turn, an Article object acts as a
factory for Author objects, a modified dictionary built from PubMed affiliation strings.

Any bugs in the the script should be reported to matthew.emery@stemcell.com as soon as possible.

At runtime, type your stemcell username and then follow the prompts. {Could be more detailed)"""

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
    """A Batch object contains a list of Article objects and field names that are combined when writing a CSV.

    A batch is the superclass of both ZoteroEntry and Connexon objects. The methods and attributes that they have in
    common can be found here. While it is poossible to inistaniate a Batch object, it is preferable to use its
    subclass's instead.

    Attributes:
        stem_email: A valid stemcell.com email. The NIH requires the use of a valid email adress to access it's
        databases.
        articles: List that contains Article objects
        pubmed_xml: A BeautifulSoup XML object that containing the results of a PubMed eSearch.
        info: A dictionary containing any Batch-level information that the user wishes to be written in a CSV
        (i.e. Newsletter Issue or Google Scholar Search Term)
        field_names:  A tuple containing headers the user wishes to see in the outgoing CSV in the order they
        will be written
    """

    def __init__(self, stem_email, info=None, field_names=None):
        self.articles = []
        self.pubmed_xml = None
        self.stem_email = stem_email
        if info:
            self.info = info
        else:
            self.info = {'Object Type': 'Batch Test'}
        if field_names:
            self.field_names = field_names
        else:
            self.field_names = ('Article Title', 'PMID', 'Last Name', 'First Name', 'Email', 'Company', 'Department')

    def add_article(self, article):
        """Adds an article to self.articles

        Adds an Article object to self.articles. Will raise AssertionError is a non-Article object is added. Also
        prints the article's title for diagnostic pruposes.

        Arguments:
            article: An Article object

        Returns:
            None. self.articles is appended with an Article object

        Raises:
            AssertionError: Only Articles go in self.articles
        """
        assert isinstance(article, Article), 'Only Articles go in self.articles'
        print article
        self.articles.append(article)

    def create_pubmed_xml(self):
        """Returns a BeautifulSoup object from a list of Pubmed IDs.

        Creates BeautifulSoup XML object from a list of Pubmed IDs. Note that in order to minimize the stress on
        NIH servers and increase speed it is assumed that with function will be called only once in the
        construction of a Batch object. If a previous eSearch for an article's PMID was unsuccessful, there will be no
        attempt to fetch that article and the user will have to do a manual search.

        Returns:
            None. The already initialized self.pubmed_xml is replaced with a BeautifulSoup object
        """
        Entrez.email = self.stem_email
        queries = [article.get_info('PMID') for article in self.articles if article.in_info('PMID')]
        handle = Entrez.efetch(db='pubmed', id=queries, retmode='xml')
        self.pubmed_xml = BeautifulSoup(handle.read())

    def parse_pubmed_soup(self):
        """Adds a BeautifulSoup Tag object to an Article object by finding it's PMID from within self.pubmed_xml.

        Requires that self.pubmed_xml exist (i.e. self.create_pubmed_xml must have already been executed. Only article
        objects in self.articles with valid PMIDs are searched. All PubMed information about the article can be found
        within its Tag. Although Tag can be found within the Article objects info dictionary and therefore
        could be added to a CSV, it is not encouraged as the Tag would be quite large.

        Returns:
            None. The each Article object with a PMID will have an attached Tag

        Raises:
            AssertionError: 'Can't parse a PubMed soup that isn\'t there. Try self.create_pubmed_xml'
        """
        if self.pubmed_xml:
            for article in self.articles:
                if article.in_info('PMID'):
                    child = self.pubmed_xml.find('pmid', text=re.compile(article.get_info('PMID')))
                    article.update_info_dict('Tag', child.find_parent('pubmedarticle'))
        else:
            raise AssertionError('Can\'t parse a PubMed soup that isn\'t there. Try self.create_pubmed_xml')

    def write_csv(self, csv_file):
        """Writes a CSV file with by combining the info dictionaries of the Batch, Articles and Author objects.

        The method leverages the powerful DictWriter object to build a highly customizable CSV out of the info
        dictionaries of each class. The appearance and order of CSV is dictated by self.field_names, which can be
        defined at initialization. DictWriter will leave empty any information it does not know.

        Arguments:
            csv_file: An OPEN .csv file in 'wb/ab' mode. Binary mode is critical to this method's success.

        Raises:
            AssertionError: "Open CSV file in proper mode! ('wb'/'ab')"
        """
        assert csv_file.mode == 'wb' or csv_file == 'ab', "Open CSV file in proper mode! ('wb'/'ab')"
        print 'Writing CSV.'
        csv_writer = csv.DictWriter(csv_file, extrasaction='ignore', fieldnames=self.field_names)
        csv_writer.writeheader()
        for article in self.articles:
            for author in article.authors:
                full_dict = dict(author.info.items() + article.info.items() + self.info.items())
                csv_writer.writerow(full_dict)

        csv_file.close()


class ZoteroEntry(Batch):
    """Subclass of Batch object that populates self.articles from a Zotero CSV.

    Zotero is a free citation manager that, among other things, allows the user to scrape citation information from
    Google Scholar. Zotero does an excllent job of this, so it would be unproductive to code a seperate, inferior
    version for our purposes. However, there is still key information missing from this process, so this class takes
    the information already known and combines it with PubMed queries to build a more complete lead entry CSV file.

    Attributes:
        zotero_csv: An CSV file exported from Zotero
        stem_email: A valid stemcell.com email. The NIH requires the use of a valid email adress to access it's
        databases.
        articles: List that contains Article objects
        pubmed_xml: A BeautifulSoup XML object that containing the results of a PubMed eSearch.
        info: A dictionary containing any Batch-level information that the user wishes to be written in a CSV
        (i.e. Google Scholar Search Term)
        field_names: A tuple containing headers the user wishes to see in the outgoing CSV in the order they
        will be written.
    """
    def __init__(self, zotero_csv, stem_email, info=None, field_names=None):
        Batch.__init__(self, stem_email)
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
        """Reads self.zotero_csv and populates the self.articles list.

        Opens self.zotero_csv as a csv DictReader object. The reader gies through each row, appending Article Title
        and URL information. Note that gathering URL information here allows greatly increases the speed of the script,
        since collection and validation of DOIs is a significant bottleneck. At this point, the method will attempt
        to find the PMID based on article name. If this unsuccessful, function builds article.authors from the CSV.
        Finally, the article is added to the self.articles list.

        Returns:
            None. Method add partially constructed articles into self.articles.
        """
        with self.zotero_csv as zotero:
            reader = csv.DictReader(zotero)
            for row in reader:
                if row['Item Type'] == 'journalArticle':
                    article = Article(info={'Article Title': row['Title'],
                                            'Publication Link': row['Url']})
                    article.update_info_dict('PMID', article.lookup_up_pmid(article.info['Article Title'],
                                                                            self.stem_email))
                    if not article.in_info('PMID'):
                        article.update_info_dict('Authors', row['Author'].split('; '))
                    self.add_article(article)

    def construct_articles(self):
        """Further constructs articles from self.read_csv after self.pubmed_xml is constructed."""
        for article in self.articles:
            if article.in_info('Tag'):
                article.update_info_dict('Publication Date', article.find_date())
                article.update_info_dict('Abstract', article.find_abstract())
                article.update_info_dict('Authors', article.find_authors())


class Newsletter(Batch):
    """Subclass of Batch object that populates self.articles from a valid Connexon URL.

    Opens the provided Connexon URL and scrapes article title and URL information from the issue. The information
    is then processed in the same way as a Batch object.

    Attributes:
        url: A valid Connexon url
        stem_email: A valid stemcell.com email. The NIH requires the use of a valid email adress to access it's
        databases.
        soup = A BeautifulSoup object created from self.url.
        articles: List that contains Article objects.
        pubmed_xml: A BeautifulSoup XML object that containing the results of a PubMed eSearch.
        info: A dictionary containing any Batch-level information that the user wishes to be written in a CSV
        (i.e. Newsletter Issue or Google Scholar Search Term)
        field_names:  A tuple containing headers the user wishes to see in the outgoing CSV in the order they
        will be written.
    """
    def __init__(self, url, stem_email, info=None, field_names=None):
        Batch.__init__(self, stem_email)
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
        """Given a URL will return a BeautifulSoup object."""
        header = {'User-Agent': "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/534.30 (KHTML, like Gecko) Ubuntu/11.04 "
                                "Chromium/12.0.742.112 Chrome/12.0.742.112 Safari/534.30"}
        request = urllib2.Request(self.url, headers=header)
        return BeautifulSoup(urllib2.urlopen(request, timeout=20))

    def parse_connexon(self):
        """Given a Connexon BeautifulSoup object, adds Article objects to self.articles.

        Parses article information from self.soup and constructs Article objects from it. This method will scrape
        Article Title and Publication Link from self.soup and attempt a PMID lookup.

        Returns:
            None. Appends Article object to self.articles.
        """
        pubs = self.soup.find_all(_find_comment)
        for pub in pubs:
            article = Article(info={'Article Title': pub.text.lstrip('\n')})
            article.update_info_dict('PMID', article.lookup_up_pmid(pub.text.lstrip('\n'), self.stem_email))
            article.update_info_dict('Publication Link', pub.find_previous('a').get('href'))
            self.add_article(article)

    def find_specific_lead_source(self):
        """Returns the name of the specific lead source."""
        title = self.soup.find('title').text
        return title.split(' - ')[1]

    def construct_articles(self):
        """Adds Publication Date and Author Information to each article based on self.pubmed_xml."""
        for article in self.articles:
            if article.in_info('Tag'):
                article.update_info_dict('Publication Date', article.find_date())
                article.update_info_dict('Authors', article.find_authors())


class Article(object):
    """An object containing a dictionary for writable information and a list of Author objects.

    The Article object has two respionsibilities. The first is to populate itself with CSV writable information for the
    Batch write_CSV method. The second is to act as a factory for Author objects.

    Arguments:
        info:
        authors:
    """
    def __init__(self, info=None):
        if info:
            self.info = info
        else:
            self.info = {}
        self.authors = []

    def in_info(self, query):
        """Returns boolean for whether query is in self.info."""
        try:
            if self.info[query]:
                return query in self.info
        except KeyError:
            return False

    def get_info(self, query):
        """Returns query's key value from self.info. Returns None if query not in info."""
        if query in self.info:
            return self.info[query]
        else:
            return None

    def lookup_up_pmid(self, pub_title, stem_email, translated=False):
        """Returns PMID for a search of the publication title.

        Querys the NIH Entrez database with publication title and returns the PMID if successful. If unsuccessful,
        this method will translate the publication title into

        Arguments:
            pub_title:
            stem_email:
            translated:
        """
        Entrez.email = stem_email
        handle = Entrez.esearch(db='pubmed', term=pub_title, retmax=10, sort='relevance')
        try:
            return Entrez.read(handle)['IdList'][0]
        except IndexError:
            if not translated:
                pub_title = self.translate_british(pub_title)
                return self.lookup_up_pmid(pub_title, stem_email, translated=True)
            else:
                print 'Could not find PMID: {}'.format(pub_title)

    @staticmethod
    def translate_british(publication_title):
        """Translates commonly Americanized words back into their British counterparts."""
        brit_dict = {'Leukemia': 'Leukaemia',
                     'Tumor': 'Tumour',
                     'Signaling': 'Signalling',
                     'α': 'alpha'}
        for brit in brit_dict.items():
            publication_title = publication_title.replace(brit[0], brit[1])
        return publication_title

    def update_info_dict(self, key, value):
        """Simple key-value update for self.info."""
        self.info[key] = value

    def get_info_items(self):
        """Returns all of self.info in the form of a tuple of two-tuples"""
        return self.info.items()

    def find_title(self):
        """Returns a string derived from the Article Tag it's title"""
        return self.info['Tag'].articletitle.text.strip().strip('.')

    def find_date(self):
        potential_tags = [self.info['Tag'].find('pubmedpubdate', {'pubstatus': 'aheadofprint'}),
                          self.info['Tag'].pubdate,
                          self.info['Tag'].find(pubstatus='medline')]
        potential_tags = [pot_tag for pot_tag in potential_tags if pot_tag]
        for tag in potential_tags:
            year, month, day = self.return_date_from_tag(tag)
            if year and month and day:
                return self.output_date(year, month, day)
        else:
            print 'Could not find date'
            return None


    @staticmethod
    def output_date(year, month, day):
        if month.isdigit():
            pubdate = time.strptime('{}{}{}'.format(month, day, year), '%m%d%Y')
        else:
            pubdate = time.strptime('{}{}{}'.format(month, day, year), '%b%d%Y')
        return time.strftime('%m/%d/%Y', pubdate)

    @staticmethod
    def return_date_from_tag(find_tag):
        try:
            year = find_tag.find('year').text.strip()
            month = find_tag.find('month').text.strip()
            day = find_tag.find('day').text.strip()
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


def make_zotero_entry(stem_email):
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
    return ZoteroEntry(source, stem_email, info=info_dict)


if __name__ == '__main__':
    EMAIL = raw_input('Input valid Stemcell email address: ')
    assert re.match(r'[\S]+@stemcell.com', EMAIL), 'Script requires valid Stemcell email address.'
    PROMPT = raw_input('Press "Z" for Zotero. Press "C" for Connexon. ')
    if PROMPT.upper() == 'Z':
        batch = make_zotero_entry(EMAIL)
        batch.read_csv()
    elif PROMPT.upper() == 'C':
        chosen_url = url_wrapper()
        batch = Newsletter(chosen_url, EMAIL)
    else:
        raise AssertionError('Invalid choice, try again')
    batch.create_pubmed_xml()
    batch.parse_pubmed_soup()
    batch.construct_articles()
    batch.write_csv(open('BatchOutput.csv', 'wb'))

    # TODO: Catch Press Release first titles
    # TODO: Search for Salesforce IDs?