__author__ = 'memery'

from Bio import Entrez, Medline
import string
import csv
import re


def fetch_ids(query):
    Entrez.email = "matthew.emery@stemcell.com"
    handle = Entrez.esearch(db='pubmed', term=query, retmax=50, sort='relevance')
    record = Entrez.read(handle)
    idlist = record["IdList"]
    handle = Entrez.efetch(db='pubmed', id=idlist, rettype='medline', retmode='text')
    records = list(Medline.parse(handle))
    outfile = open('test.txt', 'w')
    for record_source in records:
        source = clean_source(record_source)
        outfile.write('{}: {}\n'.format(source, lookup_if(source)))
    for record in records:
        outfile.write('>>-->>Publication IF>>\n')
        outfile.write(lookup_if(clean_source(record)) + '\n')
        outfile.write('>>-->>Publication Title>>\n')
        outfile.write(clean_title(record.get("TI", "?")) + '\n')
        outfile.write('>>-->>Publication Content>>\n')
        outfile.write(record.get('AB', '?') + '\n')
        outfile.write('>>-->>Publication Journal>>\n')
        outfile.write('{}\n'.format(clean_source(record)))
        outfile.write('>>-->>Publication Link Names>>\n')
        outfile.write('Abstract\n')
        outfile.write('>>-->>Publication Link URLs>>\n')
        outfile.write(clean_doi(record))


def clean_title(text):

    text = text.title()
    ignore_list = ['in', 'a', 'the', 'of', 'between', 'against', 'and', 'for', 'on', 'with', 'or', 'without', 'by',
                   'into', 'to', 'through', 'at', 'from', 'following', 'an', 'as', 'after', 'before', 'toward',
                   'towards', 'that', 'versus', 'via', 'from', 'during', 'miR', 'but', 'followed', 'siRNA', 'vs']
    for word in ignore_list:
        if re.search(r'\b{}\b'.format(word.title()), text[1:]):
            text = re.sub(r'\b{}\b'.format(word.title()), r'{}'.format(word), text)
    text = emphasize_latin(text)
    # text = capitalize_phase(text)

    return text.rstrip('.')


def clean_source(record):
    source = record.get('SO', '?')
    source = source.split('.')[0]
    return str(source)


def clean_doi(record):
    record_list = record.get('AID', '?')
    for n in record_list:
        if 'doi' in n:
            doi = n
            break
    try:
        doi = doi.split(' ')[0]
        return 'http://dx.doi.org/{}\n'.format(doi)
    except UnboundLocalError:
        return 'No DOI!\n'


def clean_abstract(record):
    pass


def lookup_if(source):
    with open(r'S:\CXN - Connexon\Connexogen\impactfactors.csv') as fin:
        csvin = csv.reader(fin)
        lookup = {row[0]: row for row in csvin}
    try:
        return str(lookup[source][2])
    except:
        return 'Not in database'


def emphasize_latin(text, title=True):
    """Add em tags to latin phrases such as 'in vivo'"""
    latin_list = ['in vivo', 'ex vivo', 'in vitro', 'in silico', 'de novo']
    for word in latin_list:
        if re.search(r'\b{}\b'.format(word.title()), text, flags=re.IGNORECASE):
            if title:
                return re.sub(word, r'<em>{}</em>'.format(word.title()), text, flags=re.IGNORECASE)
            else:
                return re.sub(word, r'<em>{}</em>'.format(word), text, flags=re.IGNORECASE)
    return text


# def capitalize_phase(text):
#     match = re.search(r'[p|P]hase [1-3]', text)
#     if match:
#         text_split = text.lower.split(' ')
#         for word in range(len(text_split)):
#             if word == 'phase':
#                 num = text_split(word + 1)
#         if 'i' in num:
#             num = len(num)
#         return re.sub(r'[p|P]hase [1-3|i*]', 'Phase {}'.format(num)), text)
#     return text


if __name__ == '__main__':
    prompt = raw_input('Pubmed query? ')
    fetch_ids(prompt)
    print 'Finished'

# TODO: Alzheimer'S
# TODO: Add EPDAT prompt
# TODO: Open test.txt