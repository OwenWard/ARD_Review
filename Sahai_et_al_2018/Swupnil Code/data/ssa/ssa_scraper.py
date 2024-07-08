import urllib
from bs4 import BeautifulSoup

#write table to csv file
f = open('ssa_names.csv','w')

url = 'https://www.ssa.gov/cgi-bin/popularnames.cgi'

#iterate over all years
for year in range(1900, 2015):
    post_data = {
        'year': str(year),
        'top': str(500),
        'number': 'p'
    }

    post_encode = urllib.urlencode(post_data)
    post_encode = post_encode.encode('UTF-8')
    page = urllib.urlopen(url, post_encode).read().decode('UTF-8', 'ignore')
        
    soup = BeautifulSoup(page)
    table = soup.find("table",{"width":"72%"})

    if year % 10 == 0 :
        print year
    elif year % 5 == 0 :
        print '.'

    #iterate over each row
    for row in table.findAll("tr"):
        cells = row.findAll("td")
        if len(cells) == 5:
            rank = cells[0].find(text=True)
            male = cells[1].find(text=True)
            maleP = cells[2].find(text=True).replace("%","")
            female = cells[3].find(text=True)
            femaleP = cells[4].find(text=True).replace("%","")
            write_to_file = str(year) + "," + rank + "," + male + "," + maleP + "," + female + "," + femaleP + "\n"
            f.write(write_to_file)

f.close()
