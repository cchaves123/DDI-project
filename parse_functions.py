from datetime import datetime
import numpy as np
import re



def parse_agent(c):
    id_c = re.findall(r'<id>(.*?)<\/id>', c)[0] if len(re.findall(r'<id>(.*?)<\/id>', c))>0 else ""
    identifier= re.findall(r'<identifier>(.*?)<\/identifier>', c)[0] if len(re.findall(r'<identifier>(.*?)<\/identifier>', c))>0 else ""
    poly_id = re.findall(r'<polypeptide id=\"(.*?)\"', c)[0] if len(re.findall(r'<polypeptide id=\"(.*?)\"', c))>0 else ""
    synonyms = re.findall(r'<synonym>(.*?)<\/synonym>', c)
    if len(synonyms) == 0: synonyms =[]
    name = re.findall(r'<name>(.*?)<\/name>', c)[0] if len(re.findall(r'<name>(.*?)<\/name>', c))>0 else ""
    action = re.findall(r'<action>(.*?)<\/action>', c)
    if len(action) == 0: action =[]
    pfams = re.findall(r'<identifier>(.*?)<\/identifier>', c)
    if len(pfams) == 0: pfams =[]
    go_classifiers = re.findall(r'<go-classifier>(.*?)<\/go-classifier>', c)
    if len(go_classifiers) == 0: go_classifiers =[]
    cellular_location = re.findall(r"<cellular-location>(.*?)<\/cellular-location>", c)
    if len(cellular_location) == 0: cellular_location =[]
    return(poly_id, {'name':name, 'id': id_c,
            'identifier': identifier,
            'cellular_location' : cellular_location,
            'action': action,
            'pfam': pfams,
            'synonyms': synonyms,
            'go-classifier': [[re.findall(r"<category>(.*?)<\/category>", p)[0], re.findall(r"<description>(.*?)<\/description>", p)[0]]
                    for p in go_classifiers]
                        })

def parse_classification(classification):
    pattern = r"<direct-parent>(.*?)<\/direct-parent>"
    direct_parent=re.findall(pattern, classification)
    if len(direct_parent) == 0:
        res = {'direct-parent' : 'None'}
    else:
        res = {'direct-parent' : direct_parent[0]}

    pattern = r"<kingdom>(.*?)<\/kingdom>"
    kingdom = re.findall(pattern, classification)
    if len(kingdom) == 0:
        res['kingdom'] = 'None'
    else:
        res['kingdom'] = kingdom[0]

    pattern = r"<superclass>(.*?)<\/superclass>"
    superclass = re.findall(pattern, classification)
    if len(superclass) == 0:
        res['superclass'] = 'None'
    else:
        res['superclass'] = superclass[0]

    pattern=r"<class>(.*?)<\/class>"
    Class=re.findall(pattern, classification)
    if len(Class)==0:
        res['class'] = 'None'
    else:
        res['class'] = Class[0]

    pattern=r"<subclass>(.*?)<\/subclass>"
    subclass=re.findall(pattern, classification)
    if len(subclass)==0:
        res['subclass'] = 'None'
    else:
        res['subclass'] = subclass[0]
    return(res)

def parse_pathways(pathlist):
    out = {}
    if len(pathlist) > 0:
        for j in range(len(pathlist)):
            index = 'pathways'+str(j+1)
            #pattern = r"<pathway>(.*)<ds>"
            path = pathlist[j]
            smpdb_id = re.findall(r"<smpdb-id>(.*)</smpdb-id>", path)[0] if len(re.findall(r"<smpdb-id>(.*)</smpdb-id>", path))>0 else ""
            namep = re.findall(r"<name>(.*?)<\/name>", str(path).replace('.', '. '))[0]
            category = re.findall(r"<category>(.*)</category>", str(path))[0] if len(re.findall(r"<category>(.*)</category>", path))>0 else ""
            pathways_details = {'smpdb_id':smpdb_id,
                                'name_pathway': namep,
                                'category_pathway' : category}

            pattern=r"<drug>(.*?) </drug>"
            druglist=re.findall(pattern, pathlist[j])
            if len(druglist)>0:
                drugs_in_pathway  = [re.findall(r"<drugbank-id>(.*)<\/drugbank-id>", d)[0] for d in druglist]
            pathways_details['drugs_in_pathway'] = drugs_in_pathway

            pattern="<enzymes>(.*)</enzymes>"
            enzyme_list = re.findall(pattern, pathlist[j])
            enzymes_in_pathway = []
            if len(enzyme_list)>0:
                enzymes_in_pathway = re.findall(r"<uniprot-id>(.*?)</uniprot-id>", enzyme_list[0])
            pathways_details['enzymes_in_pathway'] = enzymes_in_pathway
            out[index] = pathways_details
    return(out)

def parse_ddi(newlist):
    key = 1
    out = {}
    if len(newlist) > 0:
        for j in range(len(newlist)):
            index = 'drug_interaction' + str(key)
            pattern = r"<drugbank-id>(.*)<\/drugbank-id>"
            match = re.findall(pattern,newlist[j])
            drug_bank_id = match[0]
            pattern1 = r"<name>(.*)<\/name>"
            name1 = re.findall(pattern1,newlist[j])[0]
            pattern2 = r"<description>(.*)</description>"
            description1 = re.findall(pattern2,newlist[j])[0]
            out[drug_bank_id] = {'name':name1,
                                  'description': description1}
    return(out)


def parse_categories(x):
    cat_code  = re.findall(r"<categories>(.*)</categories>", x.lower())
    if len(cat_code)>0:
        all_cats = [[e[0], e[2]] for e in re.findall(r'<category>[ ]+<category>(.*?)<\/category>[ ]+(<mesh-id>|)(.*?)(<\/mesh-id>|<mesh-id\/>)[ ]+ <\/category>', cat_code[0])]
    else:
        all_cats = []
    return(all_cats)



def parse_atc_codes(x):
    atc_code  = re.findall(r"<atc-codes>(.*)</atc-codes>", x.lower())
    if len(atc_code)>0:
        code = re.findall(r'<atc-code code=\"(.*?)\">', atc_code[0])
        all_codes = re.findall(r'<level code=\"(.*?)\">(.*?)</level>', atc_code[0])
    else:
        code = ''
        if len(re.findall(r"</atc-codes>", x)):
            all_codes = []
        else:
            print("Drug has no atc code")
            all_codes = []
    return(code, all_codes)

def find_eos(x):
    usetxt = x
    #print(newtxt)
    drugbank_id = re.findall(r'<drugbank-id\ primary=\"true\">(.*?)<\/drugbank-id>',
                              x)[0]

    name = re.findall(r"<name>(.*?)<\/name>", usetxt)[0]
    description_drug = re.findall(r"<description>(.*?)<\/description>", usetxt)
    if len(description_drug) == 0:
        description_drug = ""
    else:
        description_drug = description_drug[0]

    indication = re.findall(r"<indication>(.*)<\/indication>", usetxt)
    if len(indication) == 0:
        indication = ""
    else:
        indication = indication[0]

    pattern = pattern = r"<pharmacodynamics>(.*)<\/pharmacodynamics>"
    pharmacodynamics = re.findall(pattern, usetxt)
    if len(pharmacodynamics) == 0:
        pharmacodynamics = ""
    else:
         pharmacodynamics = re.findall(pattern, usetxt)[0]

    pattern = r"<mechanism-of-action>(.*)<\/mechanism-of-action>"
    mechanism_of_action =  re.findall(pattern, usetxt)
    if len(mechanism_of_action) == 0:
        mechanism_of_action = ""
    else:
        mechanism_of_action = re.findall(pattern, usetxt)[0]

    toxicity = re.findall(r"<toxicity>(.*)<\/toxicity>", usetxt)
    if len(toxicity) == 0:
        toxicity = ""
    else:
        toxicity = toxicity[0]

    pattern = r"<metabolism>(.*)<\/metabolism>"
    metabolism = re.findall(pattern, usetxt)
    if len(metabolism) == 0:
        metabolism = ""
    else:
         metabolism = re.findall(pattern, usetxt)[0]

    absorption = re.findall(r"<absorption>(.*)<\/absorption>", usetxt)
    if len(absorption) == 0:
        absorption = ""
    else:
        absorption = absorption[0]

    half_life= re.findall(r"<half-life>(.*)<\/half-life>", usetxt)
    if len( half_life) == 0:
         half_life = ""
    else:
         half_life= half_life[0]

    affected_organism = re.findall(r"<affected-organism>(.*?)<\/affected-organism>",
                                     usetxt)

    group = re.findall(r"<group>(.*?)<\/group>",
                                     usetxt)[0]

    pattern=r"<route-of-elimination>(.*)</route-of-elimination>"
    route_of_elimination= re.findall(pattern, usetxt)
    if len(route_of_elimination) == 0:
        route_of_elimination = ""
    else:
        route_of_elimination = re.findall(pattern, usetxt)[0]


    pattern = r"<started-marketing-on>(.*?)<\/started-marketing-on>"
    market_start = re.findall(pattern, usetxt)
    if len(market_start) > 0:
        market_start = [datetime.strptime(m, '%Y-%m-%d').date().year
                for m in market_start]
        market_start = np.min(market_start)
    else:
        market_start = np.nan

    pattern = r"<ended-marketing-on>(.*?)</ended-marketing-on>"
    market_end = re.findall(pattern, usetxt)
    if len(market_end) > 0:
        market_end = [datetime.strptime(m, '%Y-%m-%d').date().year
                    for m in market_end]
        market_end = np.max(market_end)
    else:
        market_end = 2050


    pattern=r"<volume-of-distribution>(.*)</volume-of-distribution>"
    volume_of_distribution= re.findall(pattern, usetxt)
    if len( volume_of_distribution) == 0:
         volume_of_distribution = ""
    else:
        volume_of_distribution = re.findall(pattern, usetxt)[0]


    classification= re.findall(r"<classification>(.*)</classification>", usetxt)
    if len(classification) == 0:
         classification = {}
    else:
        classification = classification[0]
        classification = parse_classification(classification)

    #go_classifiers = [[re.findall(r"<category>(.*?)<\/category>", p)[0], re.findall(r"<description>(.*?)<\/description>", p)[0]]
    #                    for p in re.findall(r"<go-classifier>(.*?)<\/go-classifier>", usetxt)]


    target_str= re.findall(r"<targets>(.*?)<\/targets>", usetxt)
    targets = {}
    if len(target_str)>0:
        target= re.findall(r"<target>(.*?)<\/target>", target_str[0])
        if len(target) == 0:
             target = ""
        else:
            for c in target:
                polypeptide_id, description = parse_agent(c)
                targets[polypeptide_id] = description

    carrier_str= re.findall(r"<carriers>(.*?)<\/carriers>", usetxt)
    carriers = {}
    if len(carrier_str)>0:
        carrier_lst = re.findall(r'<carrier(.*?)<\/carrier>', carrier_str[0])
        for c in carrier_lst:
            polypeptide_id, description = parse_agent(c)
            carriers[polypeptide_id] = description

    transporters_str= re.findall(r"<transporters>(.*?)<\/transporters>", usetxt)
    transporters = {}
    if len(transporters_str)>0:
        transporters_lst = re.findall(r'<transporter(.*?)<\/transporter>', transporters_str[0])
        for c in transporters_lst:
            polypeptide_id, description = parse_agent(c)
            transporters[polypeptide_id] = description

    enzyme_str= re.findall(r"<enzymes>(.*?)<\/enzymes>", usetxt)
    if len(enzyme_str)>0:
        enzymes =  re.findall(r"<uniprot-id>(.*?)<\/uniprot-id>", enzyme_str[0])
    else:
        enzymes = []

    pattern=r"<food-interaction>(.*?)<\/food-interaction>"
    food_interactions= re.findall(pattern, usetxt)
    if len(food_interactions) == 0:
         food_interactions = ""
    else:
        food_interactions = food_interactions
    #food_interactions = parse_food_interactions(food_interactions)
    #food_interactions = parse_food_interactions(food_interactions)

    pattern=r"<drug-interaction>(.*?)</drug-interaction>"
    drug_interactions= re.findall(pattern, usetxt)
    if len(drug_interactions) == 0:
         ddi_list = {}
    else:
        ddi_list = parse_ddi(drug_interactions)

    ###

    pattern = r"<pathways>(.*?) </pathways>"
    pathways = re.findall(pattern, usetxt)
    if len(pathways) == 0:
        pathways = ""
        pathlist=""
    else:
        pathlist = re.findall(r"<pathway>(.*?) </pathway>", pathways[0])
    pathways = parse_pathways(pathlist)


    atc_code, all_atc_codes = parse_atc_codes(usetxt)
    cat_codes  = parse_categories(usetxt)

    interest = re.findall("<external-identifiers>(.*)</external-identifiers>", usetxt)
    if len(interest) > 0:
        ident = re.findall(r'<external-identifier>[ ]+<resource>(.*?)</resource>[ ]+<identifier>(.*?)</identifier>', interest[0])
        identifiers  = {k:v for (k,v)  in ident}
    else:
        identifiers  = {}

    smile = re.findall("<property>[ ]+<kind>SMILES</kind>[ ]+<value>(.*?)</value>", usetxt)

    usetxt = list1[1]
    names = []
    names_str= re.findall(r"<products>(.*?)<\/products>", usetxt)
    if len(target_str)>0:
        in_name = re.findall(r"<product>(.*?)<\/product>", names_str[0])
        if len(in_name) == 0:
             n = ""
        else:
            for c in in_name:
                 n = re.findall(r"<name>(.*?)<\/name>", c)
                 names += [n]

    return({'drugbank_id': drugbank_id,
            'name': name,
            'description':description_drug,
            'indication': indication,
            'pharmacodynamics': pharmacodynamics,
            'toxicity': toxicity,
            'metabolism' : metabolism,
            'absorption': absorption,
            'half_life': half_life,
            'market_start': market_start,
            'affected_organism': affected_organism,
            'market_end': market_end,
            'mechanism_of_action' : mechanism_of_action,
            'route_of_elimination': route_of_elimination,
            'volume_of_distribution': volume_of_distribution,
            'classification': classification,
            'carriers': carriers,
            'targets' : targets,
            'transporters': transporters,
            'enzymes': enzymes,
            'atc_code':  atc_code,
            'identifiers' : identifiers,
            'categories': cat_codes,
            'atc_hierarchy' :  all_atc_codes,
            'food_interactions': food_interactions,
            'ddi': ddi_list,
            'smile': smile,
            'names': list(np.unique(names)),
            'pathways': pathways})
