from Functions import *
from Rewiring_functions import *

Ilist=read_database()

#create annual networks in each site
networks={}
for site in Sites:
    #print(site)
    networks[site]={}
    for year in Years:
        Imat = Ilist.xs((site,year),level=["Site_ID", "Year"]).pivot(index="Plant_gen_sp",columns="Pollinator_gen_sp",values="Frequency").fillna(0)
        networks[site][year]=Imat

#select two networks 
#hay que elegir sitio y año!
site=Sites[0] #2,1/9.3/10/6 /0.4
year=Years[4]
web1=networks[site][year]
web2=networks[site][year+1]

#quantify changes
Changes_df=rewiring_new_py(web1,web2)

web_t1=web1
web_t2=web2

fig, ax = plt.subplots(figsize=(12, 8))

order=None
overlap=None
web=(web_t1+2*web_t2)#.dropna(how="all").dropna(how="all", axis=1)
    #his introduces NaNs in those with only one year, so we have to put them back in
    #core_web=web.dropna(how="all").dropna(how="all", axis=1)
core_web=web
    
pollis1= set(web_t1.columns)
plants1= set(web_t1.index)

pollis2= set(web_t2.columns)
plants2= set(web_t2.index)

common_pollis = pollis1.intersection(pollis2)
common_plants = plants1.intersection(plants2)
  
pollis_y1=pollis1-pollis2
plants_y1=plants1-plants2
pollis_y2=pollis2-pollis1
plants_y2=plants2-plants1

for p in plants_y1:
    for a in pollis_y1:
        core_web.loc[p,a]=web_t1.loc[p,a]
    for a in common_pollis:
        core_web.loc[p,a]=web_t1.loc[p,a]
for a in pollis_y1:
    for p in common_plants:
        core_web.loc[p,a]=web_t1.loc[p,a]

for p in plants_y2:
    for a in pollis_y2:
        core_web.loc[p,a]=2*web_t2.loc[p,a]
    for a in common_pollis:
        core_web.loc[p,a]=2*web_t2.loc[p,a]
for a in pollis_y2:
    for p in common_plants:
        core_web.loc[p,a]=2*web_t2.loc[p,a]
for a in pollis_y2:
    for p in plants_y1:
        core_web.loc[p,a]=-1
for a in pollis_y1:
    for p in plants_y2:
        core_web.loc[p,a]=-1
    
core_web.replace(0, np.nan, inplace=True)
#common_pollis_ordered=
common_plants_ordered=list(web.loc[list(common_plants),list(common_pollis)].sum(axis=1).sort_values(ascending=False).index)
common_pollis_ordered=list(web.loc[list(common_plants),list(common_pollis)].sum(axis=0).sort_values(ascending=False).index)
plants_y1_ordered=list(web_t1.loc[list(plants_y1),list(pollis_y1)].sum(axis=1).sort_values(ascending=False).index)
pollis_y1_ordered=list(web_t1.loc[list(plants_y1),list(pollis_y1)].sum(axis=0).sort_values(ascending=False).index)

plants_y2_ordered=list(web_t2.loc[list(plants_y2),list(pollis_y2)].sum(axis=1).sort_values(ascending=False).index)
pollis_y2_ordered=list(web_t2.loc[list(plants_y2),list(pollis_y2)].sum(axis=0).sort_values(ascending=False).index)
ordered_pollis=common_pollis_ordered+pollis_y1_ordered+pollis_y2_ordered
ordered_plants=common_plants_ordered+plants_y1_ordered+plants_y2_ordered
    
#73d2de '#9ea5a9', "#ffd23f", '#3bceac', '#0ead69'
colors = ['#d7d7d6', '#73d2de', "#ffd23f",  '#009455']  # Example colors for -1, 1, 2, 3 respectively
# Create a ListedColormap
cmap_custom = ListedColormap(colors, name='custom_cmap')
colors2 = ['#ffffff', "#ffd23f", "#ffd23f",  "#ffd23f"]  # Example colors for -1, 1, 2, 3 respectively
cmap_custom2 = ListedColormap(colors2, name='custom_cmap2')

if (order):
    ordered_plants=order[0]
    ordered_pollis=order[1]

#fig, ax = plt.subplots()
    
if overlap is not None:
    over_ordered = overlap.reindex(index=ordered_plants,columns=ordered_pollis).replace(0, np.nan)
    web_t2_ordered=web_t2.reindex(index=ordered_plants,columns=ordered_pollis)
    #cax = ax.matshow(web_t2_ordered.loc[ordered_plants,ordered_pollis],cmap=cmap_custom2)
    cax = ax.matshow(over_ordered.loc[ordered_plants,ordered_pollis],cmap="brg",alpha=0.4)
else:
    cax = ax.matshow(core_web.loc[ordered_plants,ordered_pollis], cmap=cmap_custom, vmin=0, vmax=3)

c0='#009455'
c1='#73d2de'
c2="#ffd23f"
ax.axvline(len(common_pollis) - 0.5, color=c0, linestyle='--')  # Vertical line after list1
ax.axvline(len(common_pollis) - 0.4, color=c1, linestyle='--')  # Vertical line after list1
ax.axvline(len(common_pollis) + len(pollis_y1) - 0.5, color=c1, linestyle='--')  # Vertical line after list2
ax.axvline(len(common_pollis) + len(pollis_y1) - 0.4, color=c2, linestyle='--')  # Vertical line after list2
ax.axvline(len(common_pollis) + len(pollis_y1) + len(pollis_y2)- 0.5, color=c2, linestyle='--')  # Vertical line after list2

ax.axhline(len(common_plants) - 0.5, color=c0, linestyle='--')  # Horizontal line after list1
ax.axhline(len(common_plants) - 0.4, color=c1, linestyle='--')  # Horizontal line after list1
ax.axhline(len(common_plants) + len(plants_y1) - 0.5, color=c1, linestyle='--')  # Horizontal line after list2
ax.axhline(len(common_plants) + len(plants_y1) - 0.4, color=c2, linestyle='--')  # Horizontal line after list2
ax.axhline(len(common_plants) + len(plants_y1) + len(plants_y2) -0.5, color=c2, linestyle='--')  # Horizontal line after list2
ax.set_xticks([])
ax.set_yticks([])