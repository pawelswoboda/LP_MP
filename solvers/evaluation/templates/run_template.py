from jinja2 import Environment, FileSystemLoader
import sqlite3
import sys

conn = sqlite3.connect(sys.argv[1])
#c = conn.cursor()

lineStyles = ["red","blue","green","orange","black","yellow","brown","pink","cyan"]

env = Environment(loader = FileSystemLoader('.'))
template = env.get_template(sys.argv[2])
print template.render(conn = conn, c = conn.cursor(), instance_id = str(1), dataset_id = str(1), lineStyles = lineStyles)
