import os
import csv
import zipfile
import io
import base64
from flask import Flask, render_template, request, send_file, jsonify
from rdkit import Chem
from rdkit.Chem import Draw

app = Flask(__name__)

# Configuration
SDF_FILE = os.path.join(os.path.dirname(__file__), 'data', 'molecules.sdf')
os.makedirs("temp_images", exist_ok=True)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/process', methods=['POST'])
def process_molecules():
    # Obtenir la liste des molécules et les options
    molecules_text = request.form.get('molecules_list', '')
    use_color = request.form.get('use_color', 'true') == 'true'
    
    if not molecules_text.strip():
        return jsonify({"error": "Veuillez fournir une liste de molécules"}), 400
    
    # Traiter les molécules
    molecules_cibles = [mol.strip().lower() for mol in molecules_text.split('\n') if mol.strip()]
    
    # Traiter le fichier SDF et générer les images
    results = process_sdf(SDF_FILE, molecules_cibles, use_color)
    
    # Créer un fichier CSV pour Anki
    memory_file = io.BytesIO()
    with zipfile.ZipFile(memory_file, 'w') as zf:
        # Créer le fichier CSV
        csv_data = io.StringIO()
        csv_writer = csv.writer(csv_data)
        # Format Anki: front,back
        csv_writer.writerow(["front", "back"])
        
        for result in results:
            if 'image_data' in result:
                # Ajouter l'image au ZIP
                image_filename = f"{result['filename']}.svg"
                zf.writestr(f"media/{image_filename}", result['image_data'])
                
                # Ajouter l'entrée dans le CSV (format Anki)
                csv_writer.writerow([
                    f'<img src="{image_filename}" style="max-width:100%; max-height:300px;">',
                    result['molecule_name']
                ])
        
        # Ajouter le fichier CSV au ZIP
        zf.writestr("anki_molecules.csv", csv_data.getvalue())
    
    memory_file.seek(0)
    return send_file(
        memory_file,
        mimetype='application/zip',
        as_attachment=True,
        download_name='anki_molecules.zip'
    )

def process_sdf(sdf_path, molecules_cibles, use_color=True):
    results = []
    molecules_trouvees = set()
    
    # Charger les molécules avec strictParsing=False pour ignorer les erreurs de valence
    suppl = Chem.SDMolSupplier(sdf_path, sanitize=False, strictParsing=False)
    
    for idx, mol in enumerate(suppl):
        if mol is None:
            continue
            
        try:
            # Ne pas sanitizer la molécule pour éviter les erreurs de valence
            # Liste des propriétés à vérifier pour trouver le nom de la molécule
            properties_to_check = ['GENERIC_NAME', '_Name', 'DATABASE_NAME', 'COMMON_NAME', 'DRUGBANK_ID', 'DRUG_NAME']
            
            name = None
            for prop in properties_to_check:
                if mol.HasProp(prop):
                    prop_value = mol.GetProp(prop)
                    if prop_value.strip():
                        name = prop_value
                        break
            
            if not name:
                continue
                
            # Vérifier si la molécule correspond à l'une de nos cibles
            found_match = False
            matching_molecule = None
            
            for target in molecules_cibles:
                if target == name.lower():
                    found_match = True
                    matching_molecule = target
                    break
            
            if not found_match:
                continue
                
            molecules_trouvees.add(matching_molecule)
            
            # Nettoyage du nom pour le nom de fichier
            filename = matching_molecule.strip().replace(" ", "_").replace("/", "-")
            
            # Générer une image SVG
            try:
                # Paramètres de dessin
                drawer = Draw.MolDraw2DSVG(450, 400)  # Augmenter légèrement la largeur pour plus d'espace
                
                # Personnalisation des options de dessin
                opts = drawer.drawOptions()
                if not use_color:
                    opts.useBWAtomPalette()
                
                # Réduire la taille de la police
                opts.setFontSize(0.8)  # Réduction de la taille de police
                
                # Ajuster l'espacement des liaisons
                opts.bondLineWidth = 1.8  # Épaisseur des liaisons optimisée
                opts.multipleBondOffset = 0.18  # Espacement entre les liaisons multiples
                opts.padding = 0.15  # Ajouter un peu de padding
                
                # Dessiner la molécule
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText()
                
                results.append({
                    'molecule_name': name,
                    'filename': filename,
                    'image_data': svg
                })
            except Exception as e:
                print(f"Impossible de dessiner la molécule {matching_molecule}: {e}")
                
        except Exception as e:
            print(f"Erreur avec une molécule à l'index {idx}: {e}")
            continue
    
    # Ajouter des informations sur les molécules non trouvées
    not_found = list(set(molecules_cibles) - molecules_trouvees)
    
    return results

if __name__ == '__main__':
    app.run(debug=True)