import numpy as np
import pandas as pd
import plotly
import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objs as go
megan_csv = "/Users/timolucas/Documents/spirito/final.txt"
html="/Users/timolucas/Documents/spirito/spirito_comparison_abundances_heatmap.html"
pdf = "/Users/timolucas/Documents/spirito/spirito_comparison_abundances_heatmap.pdf"
threshold = 12000
h = 1200
w = 1000
font_size = 15
# xaxis=['Reactor 1 - time point 1', 'Reactor 1 time point 2','Reactor 1 time point 3','Reactor 2 time point 1','Reactor 2 - time point 2','Reactor 2 - time point 3','Reactor 3 - time point 1','Reactor 3 - time point 2','Reactor 3 - time point 3']
xaxis=['Reactor 1 - time point 1', 'Reactor 1 time point 2','Reactor 1 time point 3','Reactor 2 time point 1','Reactor 2 - time point 2','Reactor 2 - time point 3','Reactor 3 - time point 1','Reactor 3 - time point 2', 'Reactor 3 - time point 3']
data = np.genfromtxt(megan_csv, delimiter="\t")
# print(data)

data = pd.read_csv(megan_csv,delimiter='\t',header=0)



# print(filtered)
data_without_header = pd.read_csv(megan_csv,delimiter='\t')
numeric_data = data.drop("#Datasets",axis=1)
filtered = numeric_data.loc[(numeric_data.sum(axis=1) > threshold)]
otu_names = data.loc[(numeric_data.sum(axis=1) > threshold)]

print(filtered)

otu_names = list(otu_names.iloc[:,0])
for n,x in enumerate(otu_names):
    otu_names[n] = f"<i>{x}</i>"

print(otu_names)

samples = data.head(0)
samples=samples.keys().to_list()
samples=samples[1:len(samples)]
#
# print(data)

# fig = go.Figure(data=[heatmap, lines], layout=layout)
heatmap = plotly.graph_objs.Heatmap(y=otu_names,z=filtered,x=xaxis,colorscale='reds')
layout=plotly.graph_objs.Layout(xaxis=dict(side='bottom'), width=w,
                      height=h)
# lines = plotly.graph_objs.Scatter(x=["Reactor 1 - time point 3","Reactor 1  time point 3"],
#                    y=[names[8],names[0]],
#                    mode='lines',
#                    line_color='black', line_width=2.5)

fig = plotly.graph_objs.Figure(data=[heatmap], layout=layout)

fig.update_layout(font=dict(
        family="Times new roman",
        size=font_size,
        color="black"
    ))


# fig.add_trace(go.Heatmap(z=df, y=pfam_domain_list, x=list(gene_presence_dict.keys()),zmin=0,colorscale="blackbody"))

# print(list(gene_presence_dict.keys()))
#
plotly.io.orca.config.executable = '/Users/timolucas/miniconda3/bin/orca'



plotly.offline.plot(fig, filename=html, auto_open=False)

plotly.io.write_image(fig,format="pdf",file=pdf)

# fig.update_xaxes(side="top")
fig.show()

production_values = [107.30,42.43,56.98,88.06,36.81,60.05,40.01,53.98,60.95]
layout2=plotly.graph_objs.Layout(xaxis=dict(side='bottom'), width=w,
                      height=1200)

production_fig =  go.Figure(data=go.Heatmap(x=xaxis,
                    z=[production_values],y=[""],colorscale='Blues'),layout=layout2)

production_fig_annotated =  ff.create_annotated_heatmap(x=xaxis,
                    z=[production_values],y=["Caprylate total production rate (mM - C/d)"],colorscale='Blues')
production_fig_annotated.update_layout(layout2)

production_fig_annotated.update_layout(font=dict(
        family="Times new roman",
        size=font_size,
        color="black"
    ))


production_fig_annotated.show()
production_fig.show()


