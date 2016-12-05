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

lineStyles = [
    "",
    "red,mymark={star}{solid}",
    "blue,mymark={o}{solid}",
    "green,mymark={square}{solid}",
    "orange,mymark={otimes}{solid}",
    "black,mymark={triangle}{solid}",
    "yellow,mymark={diamond}{solid}",
    "brown,mymark={x}{solid}",
    "pink,mymark={+}{solid}",
    "gray,mymark={asterisk}{solid}",
    "magenta,mymark={halfcircle}{solid}"
    ]
legendStyles = [
    "",
    "red,mark=star",
    "blue,mark=o",
    "green,mark=square",
    "orange,mark=otimes",
    "black,mark=triangle",
    "yellow,mark=diamond",
    "brown,mark=x",
    "pink,mark=+",
    "gray,mark=asterisk",
    "magenta,mark=halfcircle"
    ]
    
#lineStyles = ["","red,mark=*","blue,mark=square*","green,mark=otimes*","orange,mark=triangle*","black,mark=diamond*","yellow,mark=x","brown,mark=+","pink,mark=asterisk","gray,mark=halfcircle*","magenta,mark=pointed star"] # loop index starts at 1, therefore first entry is not needed
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

if sys.argv[2] == "bound_plot.tex":
    for dataset_row in conn.cursor().execute("SELECT id,name FROM Datasets;"):
      boundPlotFileName = sys.argv[3]
      boundPlotFileName = re.sub('\.tex$', "_" + re.sub(" ", '_', dataset_row[1]) + ".tex", boundPlotFileName);
      boundPlotFile = open(boundPlotFileName,"w")
      boundPlotFile.write(template.render(conn = conn, c = conn.cursor(), dataset_id = dataset_row[0], lineStyles = lineStyles))
      boundPlotFile.close()
      call(["pdflatex", boundPlotFileName])

if sys.argv[2] == "dataset_plot.tex":
  # first we have to extend iterations such that every solver appears to have run for the same number of iterations. Do this in tmp_iterations;
  conn.cursor().execute("DROP TABLE IF EXISTS tmp_iterations;")
  conn.cursor().execute("CREATE TABLE tmp_iterations AS SELECT * FROM Iterations;")
  max_iter = conn.cursor().execute("SELECT MAX(iteration) FROM tmp_iterations;").fetchone()[0]

  for solver_instance in conn.cursor().execute("SELECT MAX(lowerBound) as lowerBound, MIN(upperBound) as upperBound, max(runtime) as runtime, max(iteration) as iteration, instance_id, solver_id FROM Iterations GROUP BY instance_id, solver_id;"):
    max_solver_iter = solver_instance[3]
    #print max_solver_iter , " < " , max_iter
    for cur_iter in range(max_solver_iter+1, max_iter+1):
      conn.cursor().execute("INSERT INTO tmp_iterations (lowerBound,upperBound,runtime,iteration,instance_id,solver_id) VALUES (?,?,?,?,?,?)", (solver_instance[0], solver_instance[1], solver_instance[2], cur_iter, solver_instance[4], solver_instance[5]))

  conn.cursor().execute("CREATE INDEX tmp_iterations_index ON tmp_iterations(solver_id, instance_id, lowerBound, upperBound, runtime, iteration);")
  for dataset_row in conn.cursor().execute("SELECT id,name FROM Datasets;"):
    datasetPlotFileName = sys.argv[3]
    datasetPlotFileName = re.sub('\.tex$', "_" + re.sub(" ", '_', dataset_row[1]) + ".tex", datasetPlotFileName);
    datasetPlotFile = open(datasetPlotFileName,"w")
    datasetPlotFile.write(template.render(conn = conn, c = conn.cursor(), dataset_id = dataset_row[0], lineStyles = lineStyles, legendStyles = legendStyles))
    datasetPlotFile.close()
    call(["pdflatex", datasetPlotFileName])
  conn.cursor().execute("DROP TABLE tmp_iterations;")

if sys.argv[2] == "dataset_plot_runtime.tex":
  for dataset_row in conn.cursor().execute("SELECT id,name FROM Datasets;"):
    datasetPlotFileName = sys.argv[3]
    datasetPlotFileName = re.sub('\.tex$', "_" + re.sub(" ", '_', dataset_row[1]) + ".tex", datasetPlotFileName);
    datasetPlotFile = open(datasetPlotFileName,"w")
    datasetPlotFile.write(template.render(conn = conn, c = conn.cursor(), dataset_id = dataset_row[0], lineStyles = lineStyles, legendStyles = legendStyles))
    datasetPlotFile.close()
    call(["pdflatex", datasetPlotFileName])

if sys.argv[2] == "dataset_convergence_runtime_plot.tex":
  for dataset_row in conn.cursor().execute("SELECT id,name FROM Datasets;"):
    datasetPlotFileName = sys.argv[3]
    datasetPlotFileName = re.sub('\.tex$', "_" + re.sub(" ", '_', dataset_row[1]) + ".tex", datasetPlotFileName);
    datasetPlotFile = open(datasetPlotFileName,"w")
    datasetPlotFile.write(template.render(conn = conn, c = conn.cursor(), dataset_id = dataset_row[0], lineStyles = lineStyles, legendStyles = legendStyles))
    datasetPlotFile.close()
    call(["pdflatex", datasetPlotFileName])

if sys.argv[2] == "runtime_plot.tex":
    for dataset_row in conn.cursor().execute("SELECT id,name FROM Datasets;"):
        for instance_row in conn.cursor().execute("SELECT id,name FROM Instances WHERE dataset_id = '" + str(dataset_row[0]) + "';"):
            runtimePlotFileName = "runtime_" + dataset_row[1] + "_" + instance_row[1]  + ".tex"
            runtimePlotFile = open(runtimePlotFileName, "w");
            runtimePlotFile.write( template.render(conn = conn, c = conn.cursor(), instance_id = str(instance_row[0]), dataset_id = dataset_row[0], lineStyles = lineStyles, legendStyles = legendStyles) )
            runtimePlotFile.close()
            call(["pdflatex", runtimePlotFileName])

if sys.argv[2] == "iteration_plot.tex":
    for dataset_row in conn.cursor().execute("SELECT id,name FROM Datasets;"):
        for instance_row in conn.cursor().execute("SELECT id,name FROM Instances WHERE dataset_id = '" + str(dataset_row[0]) + "';"):
            iterationPlotFileName = "iteration_" + dataset_row[1] + "_" + instance_row[1]  + ".tex"
            iterationPlotFile = open(iterationPlotFileName, "w");
            iterationPlotFile.write( template.render(conn = conn, c = conn.cursor(), instance_id = str(instance_row[0]), dataset_id = dataset_row[0], lineStyles = lineStyles, legendStyles = legendStyles) )
            iterationPlotFile.close()
            call(["pdflatex", iterationPlotFileName])
