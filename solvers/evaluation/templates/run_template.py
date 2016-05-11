from jinja2 import Environment, FileSystemLoader
import sqlite3
import sys

conn = sqlite3.connect(sys.argv[1])
#c = conn.cursor()

lineStyles = ["red","blue","green","orange","black","yellow","brown","pink","cyan"]
env = Environment(loader = FileSystemLoader('.'))
templateName = sys.argv[2]
template = env.get_template(templateName)

if sys.argv[2]  == "instance_table.tex":
    print template.render(conn = conn, c = conn.cursor())
if sys.argv[2] == "dataset_table.tex":
    print template.render(conn = conn, c = conn.cursor())
if sys.argv[2] == "runtime_plot.tex":
    for dataset_row in c.execute("SELECT id,name FROM Datasets;"):
        for instance_row in conn.cursor().execute("SELECT id,name FROM Instances WHERE dataset_id = '" + dataset_row[0]|string + "';"):
            print template.render(conn = conn, c = conn.cursor(), instance_id = str(instance_row[0]), dataset_id = dataset_row[0], lineStyles = lineStyles)
if sys.argv[2] == "iteratin_plot.tex":
    for dataset_row in c.execute("SELECT id,name FROM Datasets;"):
        for instance_row in conn.cursor().execute("SELECT id,name FROM Instances WHERE dataset_id = '" + dataset_row[0]|string + "';"):
            print template.render(conn = conn, c = conn.cursor(), instance_id = str(instance_row[0]), dataset_id = dataset_row[0], lineStyles = lineStyles)
