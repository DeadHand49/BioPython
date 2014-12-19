from __future__ import unicode_literals

__author__ = 'memery'

from bs4 import BeautifulSoup
import urllib2
from Bio import Entrez, Medline
import csv
import time
import re

class Newsletter(object):

    def __init__(self, url):
        self.url = url
        self.connexon_soup = self.make_soup()
        self.pub_titles = self.parse_connexon()

    def make_soup(self):
        """Given a URL will return a BeatifulSoup of that URL

        Utilizes a header to avoid 403 Errors
        """
        header = {'User-Agent': "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/534.30 (KHTML, like Gecko) Ubuntu/11.04"
                                " Chromium/12.0.742.112 Chrome/12.0.742.112 Safari/534.30"}
        request = urllib2.Request(self.url, headers=header)
        return BeautifulSoup(urllib2.urlopen(request))


    def parse_connexon(self):
        """Given a BeautifulSoup object, returns a list of publication names"""
        pubs = self.conenxon_soup.find_all(self._find_comment)
        publication_titles = []
        for pub in pubs:
            publication_titles.append(pub.text.lstrip('\n'))
        return publication_titles

    def look_up_titles(self):
        """Returns a list of PMIDs given a list of Connexon Titles"""
        pubmed_list = []
        for publication in self.pub_titles:
                pubmed_list.append(lookup_up_title(publication))

        return pubmed_list


class Article(object):

    def __init__(self):
        pass




def _find_comment(self, tag):
    """Don't use this function directly. Is used in find_comment to pull all lines of text with #PUBLICATION TITLE
    comment"""
    try:
        return tag.has_attr('face') and tag.has_attr('size') and '#PUBLICATIONS TITLE' in tag.contents[1]
    except IndexError:
        return False







def lookup_up_title(publication_title):
    Entrez.email = "matthew.emery@stemcell.com"
    handle = Entrez.esearch(db='pubmed', term=publication_title, retmax=1)
    try:
        return Entrez.read(handle)['IdList'][0]
    except IndexError:
        # return manual_pmid(publication_title)
        # uncomment for actual running, I don't understand how mock works
        return '25430711'


def manual_pmid(publication_title):
    print 'Pubmed search for ID Failed: {}'.format(publication_title)
    pmid = raw_input('Manually input PMID here: ')
    assert len(pmid) == 8, 'Malformed PMID lol'
    return pmid


def fetch_from_pubmed(pubmed_list):
    """Returns a list of Pubmed records based on a list of PMIDs"""
    Entrez.email = "matthew.emery@stemcell.com"
    handle = Entrez.efetch(db='pubmed', id=pubmed_list, retmode='xml')
    return BeautifulSoup(handle.read()).prettify()


def parse_pubmed_soup(soup):
    for article in soup('pubmedarticle'):
        authors = parse_authors(article)
        articletitle = article.find('articletitle')
        doi = article.find(eidtype="doi").strip()

def parse_authors(article):
    author_list = []
    prev_affiliation = ''
    for author in article('author'):
        email = ''
        lastname = author.find('lastname').string.strip()
        forename = author.find('forename').string.strip()
        try:
            affiliation = author.find('affiliation').string.strip()
            email_search = re.search(r'\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,4}', affiliation)
            if email_search:
                email = email_search.group(0)
            prev_affiliation = affiliation[:]
        except AttributeError:
            affiliation = prev_affiliation
        author_list.append([lastname, forename, email, affiliation])

    return author_list

def write_record(record):
    """Yields a list pertaining to single line in the CSV"""

    pub_date = get_pub_date(record)
    pub_link = clean_doi(record)
    record_list = []

    author_count = 0
    for author in record.get('FAU'):
        last_name, first_name = name_split(author)
        try:
            institute = split_institute(record.get('AD'), author_count)
            idict = parse_institute(institute)
            email = find_email(record, last_name)
        except IndexError:
            institute = 'Malformed Record'
            email = ''
        author_count += 1
        record_list.append([last_name, first_name, email, idict.get('Company'), idict.get('Department'),
                            idict.get('City'), idict.get('State'), idict.get('Country'), idict.get('Postal'),
                            'Connexon', pub_date, pub_link, record.get('TI')])

    return record_list


def write_csv(records):
    """The main function. Writes a CSV to the path directory

    Currently writes Last Name, First, Company, Publication Date, Publication Link and Publication Title

    """
    with open('leadentry.csv', 'wb') as lead_csv:
        entry_csv = csv.writer(lead_csv)
        entry_csv.writerow(['Last Name', 'First Name', 'Email', 'Company', 'Department', 'City',
                            'State', 'Country', 'Postal Code', 'Lead Source', 'Publication Date', 'Publication Link',
                            'Publication Title'])
        for record in records:
            entry_csv.writerows(write_record(record))

    lead_csv.close()


def name_split(author):
    """Returns the first and last names of an author, given a single string.

    Note: This removes the middle initials of the author"""
    last_name, first_name = author.split(', ')
    first_name = first_name.split(' ')[0]
    return last_name, first_name


def split_institute(ad, author_count):
    institutes = ad.split('. ')
    if len(institutes) == 1:
        return ad
    else:
        return institutes[author_count]


def parse_institute(institute):
    try:
        primary_institute = institute.split('; ')[0]
        department, company, city, state_and_postal, country = primary_institute.split(', ')
        state, postal = state_and_postal.split(' ')
        institute_dict = {'Department': department,
                          'Company': company,
                          'City': city,
                          'State': state,
                          'Postal': postal,
                          'Country': country.lstrip('.')}
    except ValueError:
        institute_dict = {'Department': institute,
                          'Company': '',
                          'City': '',
                          'State': '',
                          'Postal': '',
                          'Country': ''}
    return institute_dict


def regex_search(institute, mode):
    regex_dict = {'Department': r'[A-Z ]*Department[A-Z ]*|[A-Z ]*Laboratory[A-Z ]*',
                  'Company': r'[A-Z ]*University[A-Z ]*',
                  }
    query = re.search(regex_dict.get(mode), institute, flags=re.I)
    if query:
        return query.group(0)
    elif mode == 'Company':
        return institute
    else:
        return ''


def clean_doi(record):
    """Return the DOI of an article as a string"""
    record_list = record.get('AID', '?')
    for n in record_list:
        if 'doi' in n:
            doi = n
            break
    try:
        doi = doi.split(' ')[0]
        return 'http://dx.doi.org/{}'.format(doi)
    except UnboundLocalError:
        return 'No DOI!'


def get_pub_date(record):
    """Returns the publication date of a given record.

    First will attempt the electronic publication date, then will attempt the regular publication date.
    Returns a string in the form of mm/dd/yyyy"""
    try:
        epub = time.strptime(record.get('DEP'), '%Y%m%d')
    except TypeError:
        epub = time.strptime(record.get('DA'), '%Y%m%d')
    return time.strftime('%m/%d/%Y', epub)


def find_email(record, last_name):
    """Given a last name and record, will search the Institution information for an email address contain that last
    name

    Returns the email if found or a blank string if not found"""
    email = re.search(r'\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,4}', record.get('AD'), re.I)
    if not email or last_name.lower() not in email.group(0).lower():
        return ''
    if last_name.lower() in email.group(0):
        return email.group(0)


def url_wrapper():
    """Prompts user for Connexon issue URL. Will raise AssertionError if non-standard URL added.

    Returns URL"""
    url = raw_input('Paste the URL of a Connexon website here: ')
    assert '/issue/' in url, 'URL does not point to Connexon issue.'
    return url

if __name__ == '__main__':
    # test_url = url_wrapper()
    # test_url = 'http://www.mesenchymalcellnews.com/issue/volume-6-45-dec-2/'
    # test_soup = make_soup(test_url)
    # test_pub_titles = parse_connexon(test_soup)
    # test_pubmed_list = look_up_titles(test_pub_titles)
    # test_records = fetch_from_pubmed(test_pubmed_list)
    # write_csv(test_records)
    # test_soup = fetch_from_pubmed(test_pubmed_list)
    test_soup = BeautifulSoup(open('leadentry_test.xml'))
    parse_pubmed_soup(test_soup)
