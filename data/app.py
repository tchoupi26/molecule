#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import numpy as np
from rdkit.Geometry import Point3D
import sys, time, re
try:
    import genanki
except ImportError:
    print("Veuillez installer genanki : pip install genanki")
    sys.exit(1)

def sanitize_mol(mol):
    """
    Sanitize the molecule while bypassing valence check (to avoid erreurs de valence).
    """
    try:
        # On passe SANITIZE_PROPERTIES qui inclut la vérification de valence
        Chem.SanitizeMol(
            mol,
            Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
        )
    except ValueError:
        # Si ça plante encore, on l'ignore et on continue
        pass
    return mol

def get_mol_name(mol):
    """
    Tente de récupérer le nom de la molécule dans différentes propriétés.
    """
    for prop in ('name', 'Name', '_Name'):
        if mol.HasProp(prop):
            return mol.GetProp(prop)
    return None

def main():
    parser = argparse.ArgumentParser(
        description="Génère des flashcards Anki (image 2D -> nom) depuis un SDF."
    )
    parser.add_argument(
        '-i', '--input',
        default='molecules.sdf',
        help="Fichier SDF en entrée (par défaut : molecules.sdf)"
    )
    parser.add_argument(
        '-d', '--dir',
        default='flashcards_media',
        help="Répertoire pour les images (par défaut : flashcards_media/)"
    )
    parser.add_argument(
        '--names-dir',
        default='name',
        help="Répertoire contenant .txt avec DCI par ligne (par défaut : name/)"
    )
    parser.add_argument(
        '--scale',
        type=float,
        default=100.0,
        help="Pixels par unité angström (défaut : 100)"
    )
    parser.add_argument(
        '--margin',
        type=float,
        default=1.0,
        help="Marge en angström autour de la molécule (défaut : 1.0)"
    )
    args = parser.parse_args()

    os.makedirs(args.dir, exist_ok=True)
    supplier = Chem.SDMolSupplier(args.input, sanitize=False, removeHs=False)

    # pour chaque fichier de noms
    for txt in os.listdir(args.names_dir):
        if not txt.lower().endswith('.txt'):
            continue
        # lecture des DCI
        path_txt = os.path.join(args.names_dir, txt)
        with open(path_txt, encoding='utf-8') as f:
            target_names = {l.strip().lower() for l in f if l.strip()}

        flashcards = []
        for mol in supplier:
            if mol is None:
                continue
            mol = sanitize_mol(mol)
            name = None
            for prop in mol.GetPropNames():
                val = mol.GetProp(prop).strip()
                if val.lower() in target_names:
                    name = val; break
            if not name:
                continue

            # Génère les coordonnées 2D
            AllChem.Compute2DCoords(mol)
            # rotation PCA pour aligner l’axe principal horizontalement
            conf = mol.GetConformer()
            coords = np.array([[conf.GetAtomPosition(i).x,
                                conf.GetAtomPosition(i).y]
                               for i in range(mol.GetNumAtoms())])
            centroid = coords.mean(axis=0)
            c = coords - centroid
            cov = np.cov(c, rowvar=False)
            vals, vecs = np.linalg.eigh(cov)
            principal = vecs[:, np.argmax(vals)]
            angle = -np.arctan2(principal[1], principal[0])
            R = np.array([[np.cos(angle), -np.sin(angle)],
                          [np.sin(angle),  np.cos(angle)]])
            rot = c.dot(R.T) + centroid
            for i, (x, y) in enumerate(rot):
                conf.SetAtomPosition(i, Point3D(x, y, 0))
            # calcul de la taille dynamique de l’image
            xmin, ymin = rot.min(axis=0)
            xmax, ymax = rot.max(axis=0)
            w_ang = (xmax - xmin) + 2 * args.margin
            h_ang = (ymax - ymin) + 2 * args.margin
            w_px = max(1, int(w_ang * args.scale))
            h_px = max(1, int(h_ang * args.scale))

            fn = f"{name.replace(' ', '_')}.png"
            out_path = os.path.join(args.dir, fn)
            # dessin avec dimensions adaptées
            drawer = Draw.MolDraw2DCairo(w_px, h_px)
            drawer.DrawMolecule(mol); drawer.FinishDrawing()
            with open(out_path,'wb') as img:
                img.write(drawer.GetDrawingText())

            front = f'<img src="{fn}">'
            flashcards.append((front, name))

        if not flashcards:
            print(f"Aucune molécule trouvée pour {txt}.")
            continue

        # création du paquet Anki (.apkg) nommé d'après le .txt
        deck_name = os.path.splitext(txt)[0]
        deck_id = int(time.time()*1000) % (2**63)
        model = genanki.Model(
            deck_id+1, "BasicModel",
            fields=[{'name':'Front'},{'name':'Back'}],
            templates=[{
                'name':'Card 1',
                'qfmt':'{{Front}}',
                'afmt':'{{Back}}',
            }]
        )
        deck = genanki.Deck(deck_id, deck_name)
        media = []
        for front, back in flashcards:
            deck.add_note(genanki.Note(model=model, fields=[front, back]))
            m = re.search(r'src="(.+?)"', front)
            if m:
                media.append(os.path.join(args.dir, m.group(1)))

        pkg = genanki.Package(deck)
        pkg.media_files = media
        apkg = f"{deck_name}.apkg"
        pkg.write_to_file(apkg)
        print(f"✅ Généré {len(flashcards)} notes dans {apkg}")

if __name__ == "__main__":
    main()
