import os
import csv
from rdkit import Chem
from rdkit.Chem import Draw, AllChem

# Chemin vers le fichier SDF téléchargé depuis DrugBank
SDF_FILE = "molecules.sdf"

# Créer un dossier pour stocker les images
os.makedirs("images", exist_ok=True)

# Ouvrir le fichier CSV de sortie
with open("resultat_drugbank.csv", "w", newline="", encoding="utf-8") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Image", "Molécule"])

    # Charger les molécules avec sanitize=False pour éviter les problèmes de valence
    suppl = Chem.SDMolSupplier(SDF_FILE, sanitize=False)
    count = 0

    for mol in suppl:
        if mol is None:
            continue

        try:
            # Récupérer le nom du médicament (ou DCI)
            name = mol.GetProp("GENERIC_NAME") if mol.HasProp("GENERIC_NAME") else mol.GetProp("_Name")

            if not name:
                continue

            # Nettoyage du nom pour le nom de fichier
            filename = name.strip().replace(" ", "_").replace("/", "-")

            # Utiliser AllChem pour calculer les coordonnées 2D
            try:
                AllChem.Compute2DCoords(mol)
            except:
                pass

            # Générer une image SVG avec kekulize=False pour éviter certains problèmes de valence
            img_path = f"images/{filename}.svg"
            Draw.MolToFile(mol, img_path, size=(300, 300), imageType="svg", kekulize=False)

            # Écrire dans le CSV
            writer.writerow([img_path, name])
            print(f"✅ {name} -> {img_path}")
            count += 1

        except Exception as e:
            print(f"❌ Erreur avec une molécule : {e}")
            continue

print(f"\n🎉 Fini ! {count} molécules traitées.")
