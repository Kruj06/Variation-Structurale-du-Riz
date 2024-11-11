```bash
## Create a “qc” directory in the working directory in which quality control results will be stored.
# Project Riz

Ce projet porte sur l'analyse de la diversité génétique du riz après irradiation.

## Structure du projet

- **01.RawData/** : Données brutes de séquençage.
- **qc/** : Résultats du contrôle qualité.
- **metadata/** : Informations sur les échantillons.

## Exemple de commande

```bash
# Créer un dossier pour le contrôle qualité
mkdir -p qc

# Lister les fichiers dans le répertoire de données brutes
ls -lh 01.RawData/
