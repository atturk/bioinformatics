#questo script ottiene un MANE select record da NCBI e ne estrae le regioni di interazione
#è possibile specificare il gene, l'organismo e il tipo di interazione da cercare

import re
import requests
import xml.etree.ElementTree as ET


#parametri per costruire l'url
gene = 'ACE2'
organism = 'human'
interaction = 'Interaction with SARS-CoV spike glycoprotein'

#costruzione url
org = organism.replace(' ', '+')
url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term={org}[organism]+AND+{gene}[gene]+MANE+select[keyword]'

print(f'Cerco con il seguente url: {url}')

#invia richiesta
res = requests.get(url)
xml_content = res.content

#parsing xml
root = ET.fromstring(xml_content)

#cerca il tag count che servirà a controllare se ci sono risultati e quanti sono
count_tag = root.find('Count')
if count_tag is not None:
    count=count_tag.text
else:
    print('Errore: tag "Count" non trovato')

#da stringa  a int
#deve esserci un solo risultato
count = int(count)
if count == 0:
    print("Errore: nessun risultato trovato. Controllare che il gene e l'organismo siano corretti.")
elif count > 1: 
    print(f'Errore: più di un risultato trovato ({count})')
else:
    #trova l'Id che servirà per la ricerca del record
    id_element = root.find('.//Id')
    if id_element is not None:
        print(f'Un risultato trovato. ID: {id_element.text}. Procedo con la ricerca del record...')
    else:
        print('Errore: Id non trovato')

#cerca il record
id = id_element.text
url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={id}&rettype=genbank&retmode=xml'
res = requests.get(url)
genbank_content = res.content


root = ET.fromstring(genbank_content)
# pretty_xml = ET.tostring(root, encoding='unicode', method='xml')
# print(pretty_xml)

#cerca <GBFeature> e crea lista locations che conterrà le regioni di interazione
features = root.findall('.//GBFeature')
locations = []
#itera all'interno di <GBFeature>
for feature in features:
    #cerca <GBQualifier_value>
    qualifiers = feature.findall('.//GBQualifier_value')
    #itera all'interno di <GBQualifier_value>
    for qualifier in qualifiers:
        #se <GBQualifier_value> contiene la stringa di interazione, cerca la regione (GBFeature_location) la aggiunge a locations
        if interaction in qualifier.text:
            location = feature.find('.//GBFeature_location')
            locations.append(location.text)

#conta quante regioni sono state trovate
numero_regioni=len(locations)

#stampa risultati
if numero_regioni == 0:
    print('Nessuna regione compatibile con il parametro "{interaction}" trovata.')
elif numero_regioni == 1:
    print(f'{numero_regioni} regione compatibile con il parametro "{interaction}" trovata.')
else:
    print(f'{numero_regioni} regioni compatibili con il parametro "{interaction}" trovate.')

#stampa regioni trovate
print(f'\nRegione(i) trovata(e): {locations}')

#stampa sequenza nucleotidica completa:
sequence = root.find('.//GBSeq_sequence')
if sequence is not None:
    print(f' \nSequenza nucleotidica completa: \n{sequence.text.upper()}')
else:
    print('Errore: sequenza nucleotidica non trovata')

#stringa regex per estrarre le coordinate dalla lista locations
#'\d+' significa che cerca una o più cifre consecutive
reegex = re.compile(r'\d+')

#itera nelle regioni e cerca la stringa regex
for location in locations:
    #usando findall viene creata una lista con tutte le corrispondenze all'interno di un singolo elemento di location
    #un elemento di location è una stringa del tipo 'numeroInizio..numeroFine' che rappresentano inizio e fine di una regione
    match = reegex.findall(location)
    #dovrebbero esserci quindi due corrispondenze: una per l'inizio e una per la fine, controlliamo per sicurezza
    if len(match) == 2:
        #assegna inizio e fine a variabili start e end
        start = int(match[0])
        end = int(match[1])
        #stampa parametro interazione e start ed end
        print(f'\nRegione "{interaction}": {start} - {end}')
        #stampa porzione della sequenza compresa tra start ed end
        print(f'Sequenza: {sequence.text[start-1:end].upper()}')
    else:
        #se non ci sono due corrispondenze nella lista creata da findall, stampa errore
        print(f'Errore: regione non valida: {location}')

