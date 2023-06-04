import base64
from io import BytesIO
from datetime import datetime
import matplotlib.pyplot as plt


# write_html.py usage
# from write_html import write_html
# call
# write_html(input, samples, observations, pca_figs, pca_stats)

# input: str name of original vcf file
# samples: # rows of matrix (as str)
# observations: # cols of matrix (as str)
# pca_figs: list of matplot figures (as many as you want)
# pca_stats: list of top PCs (as many as you want)
# output (optional prefix for name of output file) defaults to ./output.html


def generate_info(input, samples, observations):
    '''
    Function to generate html div containing basic report information. Includes the time
    of report generation, input file name, and basic statistics
    '''
    time = str(datetime.now())
    info = '''
        <div id="Info" class="tabcontent" style="display: content">
            <table style='margin: 10% auto 10% auto'>
                <caption>Report Metadata</caption>
                <tr>
                    <td>Input</td>
                    <td>Generation Time</td>
                    <td>Samples</td>
                    <td>Observations</td>
                </tr>
                <tr>
                    <td>{input}</td>
                    <td>{time}</td>
                    <td>{samples}</td>
                    <td>{observations}</td>
                </tr>
            </table>
        </div>
    '''.format(input=input, time=time, samples=samples, observations=observations)
    return info

def generate_pca(pca_figs, pca_eigvals, pca_eigvecs):
    '''
    Function to generate html div containing plots pertaining to PCA as well as a table
    containing any relevant statistics (for example, the top 5 PCs)
    '''
    if pca_figs is None:
        return '''<div id="PCA" class="tabcontent">
                <h3 style="text-align:center">PCA</h3>
                <p> No PCA plots or statistics were generated when running EDWARD. Either an error occurred or PCA was not selected as a command-line option
                <br>Run with &lt;insert_command_here_dont_know_it_yet&gt; --pca</p>
        </div>
        '''
    all_images = ''
    for fig in pca_figs:
        tmpfile = BytesIO()
        fig.savefig(tmpfile, format='png')
        encoded = base64.b64encode(tmpfile.getvalue()).decode('utf-8')
        img = "<img src=\'data:image/png;base64,{}\'>".format(encoded)
        all_images+='\n<br>'+img
    all_vecs = ''
    for i in range(len(pca_eigvals)):
        mag = str(pca_eigvals[i])
        dir = str(pca_eigvecs[:, i])
        table_row = '''\n<tr>
                            <td>{dir}</td>
                            <td>{mag}</td>
                        </tr>'''.format(dir=dir, mag=mag)
        all_vecs += table_row


    return '''<div id="PCA" class="tabcontent">
                <h3 style="text-align:center">PCA</h3>
                <div style='float:right; margin-right: 5%'>
                    <table>
                        <caption>Top Principal Components</caption>
                        <tr>
                            <th>Direction</th>
                            <th>Magnitude</th>
                        </tr>
                        {vecs}
                    </table>
                </div>
                <div style='float:left'>
                    {imgs}
                </div>
        </div>'''.format(imgs=all_images,vecs=all_vecs)


def write_html(input,
               samples, observations,
               pca_figs=[], pca_eigvals=[], pca_eigvecs=[],
               tsne_figs=None, tsna_stats=None,
               umap_figs=None, umap_stats=None,
               output=None):
    '''
    Function to write html summary report. Expects two lists (one of figures and one of numerical statistics)
    for each of the three dimensionality reduction techniques. Unless an output prefix is specified, default to ./output.html
    '''
    # TODO ^^edit that

    # Template info div
    info_div = generate_info(input, samples, observations)
    pca_div = generate_pca(pca_figs, pca_eigvals, pca_eigvecs)

    # Template for html to insert statistics and figures into.
    # Note that {} are replaced with {{}} to be escaped (allow for use of .format)
    template = '''<html>
    <head>
        <style>
            /* Style for options menu */
            * {{ font-family: 'Franklin Gothic Medium', 'Arial Narrow', Arial, sans-serif; }}
            * {{box-sizing: border-box}}
            * {{background: #EDEAEA}}
            .tab {{
                float: left;
                border: 5px solid #ccc;
                background-color: #f1f1f1;
                width: 10%;
                height: 80%;
            }}

            /* Style for option buttons */
            .tab button {{
                display: block;
                background-color: inherit;
                color: black;
                padding: 22px 16px;
                width: 100%;
                border: none;
                outline: none;
                text-align: left;
                cursor: pointer;
                transition: 0.3s;
            }}

            /* Change background color of buttons on hover */
            .tab button:hover {{
                background-color: #ddd;
            }}

            /* Create an active/current "tab button" class */
            .tab button.active {{
                background-color: #ccc;
            }}

            /* Style the tab content */
            .tabcontent {{  
                float: left;
                padding: 0px 12px;
                border: 3px solid #ccc;
                width: 90%;
                border-left: none;
                border-bottom: none;
                height: 80%;
            }}

            img {{
                border: 2px solid#ccc;
            }}
        </style>
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>EDWARD Report</title>
    </head>

    <body>
        <style>
            div.tabcontent {{display: none}}
            table {{
                border-collapse: collapse;
                width: fit-content;
            }}

            caption {{
                background-color: #1fbbf9;
            }}

            td, th {{
                border: 1px solid #dddddd;
                text-align: center;
                padding: 8px;
                background-color: #ffffff;
            }}

        </style>

        <script>
            //Default to not showing any tabs
            tabcontent = document.getElementsByClassName("tabcontent");
            for (i = 0; i < tabcontent.length; i++) {{
                tabcontent[i].style.display = "none";
            }}

            //function to open tab on click
            function openOpt(evt, optName) {{
            var i, tabcontent, tablinks;
            tabcontent = document.getElementsByClassName("tabcontent");
            for (i = 0; i < tabcontent.length; i++) {{
                tabcontent[i].style.display = "none";
            }}
            tablinks = document.getElementsByClassName("tablinks");
            for (i = 0; i < tablinks.length; i++) {{
                tablinks[i].className = tablinks[i].className.replace(" active", "");
            }}
            document.getElementById(optName).style.display = "block";
            evt.currentTarget.className += " active";
            }}
        </script>

        <h1 style="text-align:center">EDWARD Report</h1>

        <div class="tab">
            <button class="tablinks" onclick="openOpt(event, 'Info')">Info</button>
            <button class="tablinks" onclick="openOpt(event, 'PCA')">PCA</button>
            <button class="tablinks" onclick="openOpt(event, 't-SNE')">t-SNE</button>
            <button class="tablinks" onclick="openOpt(event, 'UMAP')">UMAP</button>
        </div>
        
        {info}
        {pca}
    </body>
    </html>
    '''.format(info=info_div, pca=pca_div)

    if output is None:
        output='./output.html'
    else:
        output+='.html'
    with open(output, 'a') as file:
        file.truncate(0)
        file.write(template)
    return


## Stuff to test my html writing
input = 'fake.vcf'
samples = str(4800)
observations = str(6500)
output = './temp'
pca_figs = []
# for i in range(8):
#     fig = plt.figure()
#     plt.title("Chart #" + str(i))
#     pca_figs.append(fig)

write_html(input=input, samples=samples, observations=observations, output=output)
#################################
