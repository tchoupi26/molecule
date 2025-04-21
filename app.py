import os
import csv
from rdkit import Chem
from rdkit.Chem import Draw

# Chemin vers le fichier SDF téléchargé depuis DrugBank
SDF_FILE = "open structures.sdf"

# Créer un dossier pour stocker les images
os.makedirs("images", exist_ok=True)

# Lire la liste des molécules à traiter depuis le fichier txt
molecules_cibles = []
with open("molecules.txt", "r") as f:
    for line in f:
        molecule = line.strip()
        if molecule:
            molecules_cibles.append(molecule.lower())  # Convertir en minuscules pour comparaison insensible à la casse

print(f"🔍 Recherche de {len(molecules_cibles)} molécules spécifiques...")

# Ouvrir le fichier CSV de sortie
with open("resultat_drugbank.csv", "w", newline="", encoding="utf-8") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Image", "Molécule", "Nom dans la base"])

    # Charger les molécules avec l'option sanitize=False pour éviter les erreurs
    suppl = Chem.SDMolSupplier(SDF_FILE, sanitize=False)
    count = 0
    molecules_trouvees = set()

    for idx, mol in enumerate(suppl):
        if mol is None:
            continue

        try:
            # Essayer de sanitize la molécule, mais continuer même en cas d'échec
            try:
                Chem.SanitizeMol(mol)
            except:
                pass
            
            # Liste des propriétés à vérifier pour trouver le nom de la molécule
            properties_to_check = ['GENERIC_NAME', '_Name', 'DATABASE_NAME', 'COMMON_NAME', 'DRUGBANK_ID', 'DRUG_NAME']
            
            name = None
            # Vérifier chaque propriété
            for prop in properties_to_check:
                if mol.HasProp(prop):
                    prop_value = mol.GetProp(prop)
                    if prop_value.strip():
                        name = prop_value
                        break
            
            # Si aucun nom n'a été trouvé, passer à la molécule suivante
            if not name:
                continue

            # Vérifier si la molécule correspond EXACTEMENT à l'une de nos cibles
            found_match = False
            matching_molecule = None
            
            for target in molecules_cibles:
                # Correspondance exacte (insensible à la casse)
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
            img_path = f"images/{filename}.svg"
            try:
                Draw.MolToFile(mol, img_path, size=(300, 300), imageType="svg")
                # Écrire dans le CSV
                writer.writerow([img_path, matching_molecule, name])
                print(f"✅ Trouvé: {matching_molecule} ('{name}' dans la base) -> {img_path}")
                count += 1
            except:
                print(f"⚠️ Impossible de dessiner la molécule {matching_molecule} ('{name}')")

        except Exception as e:
            print(f"❌ Erreur avec une molécule à l'index {idx}: {e}")
            continue

    # Afficher les molécules non trouvées
    molecules_manquantes = set(molecules_cibles) - molecules_trouvees
    if molecules_manquantes:
        print("\n⚠️ Molécules non trouvées dans le fichier SDF :")
        for mol in molecules_manquantes:
            print(f"  - {mol}")

print(f"\n🎉 Fini ! {count}/{len(molecules_cibles)} molécules traitées.")
