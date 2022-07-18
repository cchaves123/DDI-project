import re


def parse_classification(classification):
    pattern = r"<direct-parent>.*</direct-parent>"
    direct_parent=re.findall(pattern, classification)
    if len(direct_parent) == 0:
        res = {'direct-parent' : 'None'}
    else:
        res = {'direct-parent' : direct_parent[0][15:-16]}

    pattern = r"<kingdom>.*</kingdom>"
    kingdom = re.findall(pattern, classification)
    if len(kingdom) == 0:
        res['kingdom'] = 'None'
    else:
        res['kingdom'] = kingdom[0][9:-10]

    pattern = r"<superclass>.*</superclass>"
    superclass = re.findall(pattern, classification)
    if len(superclass) == 0:
        res['superclass'] = 'None'
    else:
        res['superclass'] = superclass[0][12:-13]

    pattern=r"<class>.*</class>"
    Class=re.findall(pattern,classification)
    if len(Class)==0:
        res['class'] = 'None'
    else:
        res['class'] = Class[0][7:-8]

    pattern=r"<subclass>.*</subclass>"
    subclass=re.findall(pattern,classification)
    if len(subclass)==0:
        res['subclass'] = 'None'
    else:
        res['subclass'] = subclass[0][10:-11]
    return(res)


def parse_food_interactions(food_interactions, key4=1):
    foodlist=str(food_interactions).split('</food-interaction>')[:-1]
    food_int = []
    if len(foodlist)==0:
        dictionary[index_name]['food_interactions']='None'
    else:
        for j in range(len(foodlist)):
            indexf='food_interaction'+str(key4)
            food_int += [foodlist[j][23:-1]]
            key4=key4+1
    return(food_int, key4)


def parse_pathways(pathlist):
    out = {}
    if len(pathlist) > 0:
        for j in range(len(pathlist)):
            index = 'pathways'+str(j+1)
            pattern = r"<pathway>.*<drugs>"
            path = re.findall(pattern, pathlist[j])
            pattern = r"<smpdb-id>.*</smpdb-id>"
            smpdb_id = re.findall(pattern, str(path))[0][10:-11]
            #dictionary[index_name]['pathways'][index]['smpdb_id']=smpdb_id
            pattern = r"<name>.*</name>"
            namep = re.findall(pattern, str(path))[0][6:-7]
            #dictionary[index_name]['pathways'][index]['name']=namep
            pattern  =r"<category>.*</category>"
            category = re.findall(pattern, str(path))[0][10:-11]
            pathways_details = {'smpdb_id':smpdb_id,  'name of pathway': namep,
                                'category of pathway' : category}

            pattern=r"<drugs>.* </drugs>"
            drug=re.findall(pattern, pathlist[j])[0][7:-8]
            druglist=str(drug).split('</drug>')[:-1]
            key3=1
            drugs_in_pathway  = []
            for i in range(len(druglist)):
                #indexn='drugs'+str(key3)
                pattern = r"<drugbank-id>.*</drugbank-id>"
                drugbank_id = re.findall(pattern,druglist[i])[0][13:-14]
                drugs_in_pathway += [drugbank_id]
            pathways_details['drugs_in_pathway'] = drugs_in_pathway

            pattern="<enzymes>(.*)(?=</enzymes>)"
            enzymes=re.findall(pattern,pathlist[j])
            enzymes_in_pathway  = []
            if len(enzymes)>0:
                pattern="<enzymes>(.*)</enzymes>"
                enzymes_in_pathway =  [ m for i in range(len(enzyme_list)) for m in re.split('[<>]', enzyme_list[i]) if "P" in m ]
            pathways_details['enzymes_in_pathway'] = enzymes_in_pathway
            out[index] = pathways_details
    return(out)

def parse_ddi(newlist):
    key=1
    out = {}
    if len(newlist) > 0:
        for j in range(len(newlist)):
            index = 'drug_interaction' + str(key)
            pattern = r"<drugbank-id>.*</drugbank-id>"
            match = re.findall(pattern,newlist[j])
            drug_bank_id = match[0][13:-14]
            pattern1 = r"<name>.*</name>"
            name1 = re.findall(pattern1,newlist[j])
            if len(name1)==0:
                name1=''
            else:
                name1 = re.findall(pattern1,newlist[j])[0][6:-7]
            pattern2 = r"<description>.*</description>"
            description1 = re.findall(pattern2,newlist[j])[0][13:-14]
            out[name1] = description1
    return(out)


def find_eos(x):
    usetxt = x.replace('.','. ')
    #print(newtxt)
    pattern = r"<drugbank-id primary=\"true\">.*</drugbank-id>"
    drugbank_id = re.findall(pattern, usetxt)[0][28:-14]

    pattern = r"<name>.*<"
    name = re.findall(pattern, usetxt)[0][6:-1]

    pattern = r"<description>.*<"
    description = re.findall(pattern, usetxt)
    if len(description) == 0:
        description = ""
    else:
        description = re.findall(pattern, usetxt)[0][13:-1]

    pattern = r"<indication>.*<"
    indication = re.findall(pattern, usetxt)
    if len(indication) == 0:
        indication = ""
    else:
        indication = re.findall(pattern, usetxt)[0][12:-1]
    pattern = r"<pharmacodynamics>.*<"
    pharmacodynamics = re.findall(pattern, usetxt)
    if len(pharmacodynamics) == 0:
        pharmacodynamics = ""
    else:
         pharmacodynamics = re.findall(pattern, usetxt)[0][18:-1]
    pattern = r"<mechanism-of-action>.*<"
    mechanism_of_action =  re.findall(pattern, usetxt)
    if len(mechanism_of_action) == 0:
        mechanism_of_action = ""
    else:
        mechanism_of_action = re.findall(pattern, usetxt)[0][21:-1]

    pattern = r"<toxicity>.*<"
    toxicity = re.findall(pattern, usetxt)
    if len(toxicity) == 0:
        toxicity = ""
    else:
        toxicity = re.findall(pattern, usetxt)[0][10:-1]

    pattern = r"<metabolism>.*<"
    metabolism = re.findall(pattern, usetxt)
    if len(metabolism) == 0:
        toxicity = ""
    else:
         metabolism = re.findall(pattern, usetxt)[0][12:-1]

    pattern = r"<absorption>.*<"
    absorption = re.findall(pattern, usetxt)
    if len(absorption) == 0:
        absorption = ""
    else:
        absorption = re.findall(pattern, usetxt)[0][12:-1]

    pattern=r"<half-life>.*<"
    half_life= re.findall(pattern, usetxt)
    if len( half_life) == 0:
         half_life = ""
    else:
         half_life= re.findall(pattern, usetxt)[0][11:-1]
    pattern=r"<route-of-elimination>.*<"

    route_of_elimination= re.findall(pattern, usetxt)
    if len(route_of_elimination) == 0:
        route_of_elimination = ""
    else:
        route_of_elimination = re.findall(pattern, usetxt)[0][22:-1]


    pattern = r"<started-marketing-on>.*</started-marketing-on>"
    market_start = re.findall(pattern, usetxt)
    if len(market_start) > 0:
        market_start = [datetime.strptime(m[len("<started-marketing-on>"):(-len("<\started-marketing-on>"))], '%Y-%m-%d').date().year
                for m in market_start]
        market_start = np.min(market_start)
    else:
        market_start = np.nan

    pattern = r"<ended-marketing-on>.*</ended-marketing-on>"
    market_end = re.findall(pattern, usetxt)
    if len(market_end) > 0:
        market_end = [datetime.strptime(m[len("<ended-marketing-on>"):(-len("<\ended-marketing-on>"))], '%Y-%m-%d').date().year
                    for m in market_end]
        market_end = np.max(market_end)
    else:
        market_end = 2050

    newtxt = usetxt.replace('\n','.')
    pattern=r"<volume-of-distribution>.*</volume-of-distribution>"
    volume_of_distribution= re.findall(pattern, newtxt)
    if len( volume_of_distribution) == 0:
         volume_of_distribution = ""
    else:
        volume_of_distribution = re.findall(pattern, newtxt)[0][24:-25]

    pattern = r"<classification>.*</classification>"

    classification= re.findall(pattern, newtxt)
    if len(classification) == 0:
         classification = ""
    else:
        classification = re.findall(pattern, newtxt)[0][16:-17]

    pattern=r"<patents>.*</patents>"
    patents= re.findall(pattern, newtxt)
    if len(patents) == 0:
         patlist = ""
    else:
        patients = re.findall(pattern,newtxt)[0][10:-11]
        patlist = str(patents).split('</patent>')[:-1]
    pattern=r"<food-interactions>.* </food-interactions>"

    food_interactions= re.findall(pattern, newtxt)
    if len(food_interactions) == 0:
         food_interactions = ""
    else:
        food_interactions = re.findall(pattern, newtxt)[0][19:-20]
    food_interactions = parse_food_interactions(food_interactions)

    pattern=r"<drug-interactions>.*</drug-interactions>"
    drug_interactions= re.findall(pattern, newtxt)
    if len(drug_interactions) == 0:
         ddi_list = ""
    else:
        drug_interactions = re.findall(pattern, newtxt)[0][19:-20]
        ddi_list = str(drug_interactions).split('</drug-interaction>')[:-1]
        ddi_list = parse_ddi(ddi_list)

    ###

    pattern = r"<pathways>.* </pathways>"
    pathways = re.findall(pattern, newtxt)
    if len(pathways) == 0:
        pathways = ""
        pathlist=""
    else:
        pathways = re.findall(pattern, newtxt)[0][10:-11]
        pathlist = str(pathways).split('</pathway>')[:-1]
    pathways = parse_pathways(pathlist)
    return({'drugbank_id': drugbank_id,
            'name': name,
            'description':description,
            'indication': indication,
            'pharmacodynamics': pharmacodynamics,
            'toxicity': toxicity,
            'metabolism' : metabolism,
            'absorption': absorption,
            'half_life': half_life,
            'market_start': market_start,
            'market_end': market_end,
            'mechanism_of_action' : mechanism_of_action,
            'route_of_elimination': route_of_elimination,
            'volume_of_distribution': volume_of_distribution,
            'classification': classification,
            'patent_list' : patlist,
            'food_interactions': food_interactions,
            'ddi': ddi_list,
            'pathways': pathways})
