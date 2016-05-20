from jinja2 import Environment, FileSystemLoader
import sqlite3
import sys
from subprocess import call

import re

def tex_escape(text):
    """
        :param text: a plain text message
        :return: the message escaped to appear correctly in LaTeX
    """
    conv = {
        '&': r'\&',
        '%': r'\%',
        '$': r'\$',
        '#': r'\#',
        '_': r'\_',
        '{': r'\{',
        '}': r'\}',
        '~': r'\textasciitilde{}',
        '^': r'\^{}',
        '\\': r'\textbackslash{}',
        '<': r'\textless',
        '>': r'\textgreater',
    }
    regex = re.compile('|'.join(re.escape(unicode(key)) for key in sorted(conv.keys(), key = lambda item: - len(item))))
    return regex.sub(lambda match: conv[match.group()], text)


conn = sqlite3.connect(sys.argv[1])

lineStyles = ["","red","blue","green","orange","black","yellow","brown","pink","gray","magenta"] # loop index starts at 1, therefore first entry is not needed
env = Environment(loader = FileSystemLoader('.'))
env.globals['tex_escape'] = tex_escape
templateName = sys.argv[2]
template = env.get_template(templateName)

if sys.argv[2]  == "instance_table.tex":
    instanceTableFileName = sys.argv[3]
    instanceTableFile = open(instanceTableFileName,"w")
    instanceTableFile.write(template.render(conn = conn, c = conn.cursor()))
    instanceTableFile.close()
    call(["pdflatex", instanceTableFileName])

if sys.argv[2] == "dataset_table.tex":
    datasetTableFileName = sys.argv[3]
    datasetTableFile = open(datasetTableFileName,"w")
    datasetTableFile.write(template.render(conn = conn, c = conn.cursor()))
    datasetTableFile.close()
    call(["pdflatex", datasetTableFileName])

if sys.argv[2] == "runtime_plot.tex":
    for dataset_row in conn.cursor().execute("SELECT id,name FROM Datasets;"):
        for instance_row in conn.cursor().execute("SELECT id,name FROM Instances WHERE dataset_id = '" + str(dataset_row[0]) + "';"):
            runtimePlotFileName = "runtime_" + dataset_row[1] + "_" + instance_row[1]  + ".tex"
            runtimePlotFile = open(runtimePlotFileName, "w");
            runtimePlotFile.write( template.render(conn = conn, c = conn.cursor(), instance_id = str(instance_row[0]), dataset_id = dataset_row[0], lineStyles = lineStyles) )
            runtimePlotFile.close()
            call(["pdflatex", runtimePlotFileName])

if sys.argv[2] == "iteration_plot.tex":
    for dataset_row in conn.cursor().execute("SELECT id,name FROM Datasets;"):
        for instance_row in conn.cursor().execute("SELECT id,name FROM Instances WHERE dataset_id = '" + str(dataset_row[0]) + "';"):
            iterationPlotFileName = "iteration_" + dataset_row[1] + "_" + instance_row[1]  + ".tex"
            iterationPlotFile = open(iterationPlotFileName, "w");
            iterationPlotFile.write( template.render(conn = conn, c = conn.cursor(), instance_id = str(instance_row[0]), dataset_id = dataset_row[0], lineStyles = lineStyles) )
            iterationPlotFile.close()
            call(["pdflatex", iterationPlotFileName])
