

import re

def pi(annots):
    if all((X.call_rate for X in annots)):
        return(sum([X.nucl_diversity for X in annots if X.nucl_diversity != None]))



def find_SSR(st,past_ind=0):
    
    result = {}
    m = re.search(r'(.+)\1+',st) # Finds longest repeating string from start
    
    if m:
        
        i,j = m.span() # Gets the said string's index
        sub = st[i:j] # Gets the string
        ind = (sub+sub).find(sub, 1) # Finds the length of the smalest patern
        sub = sub[:ind] # Gets said smalest patern
        
        if len(sub)>1:
            if not sub in result.keys():
                result[sub] = []
            result[sub].append(st.count(sub))
            
        r = find_SSR(st[j:], past_ind+j)
        
        if r:
            for key in r:
                if key in result:
                    result[key].extend(r[key])
                else:
                    result[key] = r[key]
            return(result) # Recurcive call
        
        else:
            return(result)
    
    else:
        return(result)
    
    
    
def SSR_count(annots):
    d = {}
    for rec in annots:
        for alt in rec.ALT:
            alt = str(alt)
            d2 = find_SSR(alt)
            if len(d2):
                key = '{0}:{1}-{2}'.format(rec.CHROM, rec.start, rec.end)
                if key not in d.keys():
                    d[key] = d2
                else:
                    for key2 in d2:
                        if key2 in d[key].keys():
                            d[key][key2].extend(d2[key2])
                        else:
                            d[key][key2] = d2[key2]
    return(d)


def min_call_rate(annots):
    return(min([X.call_rate for X in annots]))