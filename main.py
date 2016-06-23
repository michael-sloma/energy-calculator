from flask import Flask
from flask import request
from flask import render_template
import stem_energy

app = Flask(__name__)
app.config['DEBUG'] = False

@app.route('/')
def stem_energy_get():
    return render_template("sequence_input.html")

@app.route('/', methods=['POST'])
def stem_energy_post():
    sequence = request.form['sequence']
    structure = request.form['structure']
    try:
        processed_text = stem_energy.energy(sequence,structure)
        return render_template("sequence_input.html", sequence=sequence, structure=structure, energy=processed_text)
    except:
        error="something went wrong; couldn't score structure. make sure sequence and structure are valid, dot-brackets are balanced, etc."
        return render_template("sequence_input.html", sequence=sequence, structure=structure, error=error)

if __name__ == '__main__':
    app.run()
