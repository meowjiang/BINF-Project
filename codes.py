import xml.etree.cElementTree as ElementTree
ElementTree.register_namespace("","http://www.drugbank.ca")
ElementTree.register_namespace("drug","http://www.drugbank.ca")
ElementTree.register_namespace("drugbank-id","http://www.drugbank.ca")

source="full_database.xml"
drugbank_ids ="drugbankIDS.txt"
context = ElementTree.iterparse(source, events=("start", "end"))
context = iter(context)

#Use this to cleanup
#event, root = context.next()


codes = dict()
for event, elem in context:
    if elem.tag == "{http://www.drugbank.ca}drug":
        for node in elem:
            if node.tag=="{http://www.drugbank.ca}drugbank-id" and "primary"in node.attrib :  
                db_id =node.text
            if node.tag=="{http://www.drugbank.ca}atc-codes":
                for code in node:
                    atc_code=code.attrib["code"]
        if event=="end":
            codes[db_id]=atc_code
            elem.clear()



file = open(drugbank_ids)
data = file.read()
lines = data.splitlines()


outfile = open("atc_codes.txt","w")
for line in lines:
    atc = codes.get(line.strip())
    outfile.write(line)
    outfile.write("\t")
    if atc==None:
        atc="ruh-roh"
    outfile.write(atc)
    outfile.write("\n")
outfile.close()

