
import re
import json



fnam = 'all_results.json'
ex = '2ovq2extract.txt'
# ex = '3fga2extract.txt'

with open(fnam) as js:
    d = json.load(js)


reslist = list()
with open(ex) as fh:
    for l in fh:
        r = l.split()
        if not r[0] in d:
            reslist.append('Protein {} not found.'.format(r[0]))
        elif not r[1] in d[r[0]]:
            reslist.append('Dimer {} not found.'.format(r[1]))
        elif not r[2] in d[r[0]][r[1]]:
            reslist.append('Mut_chain {} not found.'.format(r[2]))
        else:
            try:
                resn = int(r[3])
                all_resn = list(d[r[0]][r[1]][r[2]].keys())
                full_resns = [rr for rr in all_resn if r[3] == re.sub('[a-zA-Z]', '', rr)]
                if len(full_resns) != 1:
                    reslist.append('Mut_resnumb {} could not be mapped.'.format(r[3]))
                elif not r[4] in d[r[0]][r[1]][r[2]][full_resns[0]]:
                    reslist.append('New residue {} not found.'.format(r[4]))
                else:
                    ddG_dimer = d[r[0]][r[1]][r[2]][full_resns[0]][r[4]]
                    ddG_monomer = d[r[0]][r[2]][r[2]][full_resns[0]][r[4]]
                    if abs(ddG_dimer) > 15 or abs(ddG_monomer) > 15:
                        ddG_interface = 'NA'
                    else:
                        ddG_interface = ddG_dimer - ddG_monomer
                    reslist.append(ddG_interface)
            except:
                if not r[3] in d[r[0]][r[1]][r[2]]:
                    reslist.append('Mut_resnumb {} not found.'.format(r[3]))
                elif not r[4] in d[r[0]][r[1]][r[2]][r[3]]:
                    reslist.append('New residue {} not found.'.format(r[4]))
                else:
                    ddG_dimer = d[r[0]][r[1]][r[2]][r[3]][r[4]]
                    ddG_monomer = d[r[0]][r[2]][r[2]][r[2]][r[4]]
                    if abs(ddG_dimer) > 15 or abs(ddG_monomer) > 15:
                        ddG_interface = 'NA'
                    else:
                        ddG_interface = ddG_dimer - ddG_monomer
                    reslist.append(ddG_interface)

reslist = list(map(str, reslist))

print('\n'.join(reslist))

