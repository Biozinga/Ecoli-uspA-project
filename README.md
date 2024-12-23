# Ecoli-uspA-project (V1.0)


## Fonctionnalités

Ce programme est un **petit projet bioinformatique** permettant de :

- Vérifier la présence d'un gène d'intérêt avec un pourcentage d'identité défini.
- Vérifier la présence de boîtes consensus en amont du gène.
- Rechercher le gène *uspA* dans le génome d’*E. coli* K12.
- Générer un génome artificiel (pour comparaison statistique).
- Détecter des motifs dans la région promotrice et évaluer leur pertinence en les comparant à un génome aléatoire.


## Remarques

Le logiciel **nécessite encore des améliorations**, notamment pour :
- La **gestion et le traitement des motifs**, car beaucoup de **doublons** sont actuellement générés.
- Le fonctionnement du programme **recherche_motifs**, qui peut encore être nettement amélioré.


## Structure du projet

- `src/` : contient les fichiers sources en C :
  - `main.c`
  - `sequence_aleatoire.c`
  - `recherche_gene.c`
  - `recherche_consensus_box.c`
  - `recherche_motifs.c`
  - `assert_projet.c`
- `include/` : contient le fichier d'en-tête `projet.h`.
- `bin/` : contiendra l’exécutable final `projet_bioinfo` après compilation.
- `data/` : contient :
  - Le génome complet (`sequence_reelle.fna`).
  - La sortie des motifs détectés (`motif_retenu.txt`).
  - Le génome aléatoire (`sequence_aleatoire.fna`).
  - Le fichier du gène d'intérêt (`gene.fna`).
- `obj/` : contient les fichiers objets (.o) générés pendant la compilation.
- `README.md` : ce document.
- `Makefile` : fichier pour automatiser la compilation du projet.
- `Resultat_alignement_blastn/` : dossier servant de sauvegarde pour les résultats destinés à la publication du projet.


## Paramétrage

Le programme principal est paramétrable grâce à des variables globales définies dans **`projet.h`** (via `#define`) :

- `LONGUEUR_K_UPLET` : longueur des k-uplets analysés.
- `LONGUEUR_MIN_MOTIF` : longueur minimale d'un motif retenu.
- `X_HIT_AVANT_ENTRE` : nombre minimal d’occurrences d’un k-uplet avant extension.
- `LONGUEUR_REGION_ETUIE_MOTIFS` : longueur de la région promotrice étudiée.
- `IDENTITE_MIN` : pourcentage minimal d'identité pour détecter un gène.
- `BOITE_35` : séquence de la boîte consensus -35.
- `BOITE_10` : séquence de la boîte consensus -10.

## Compilation et exécution

1. **Se placer dans le dossier racine du projet**.
2. Exécuter la commande suivante :
   ```bash
   make
Cela générera un exécutable nommé **`projet_bioinfo`** dans le dossier **`bin/`**.

3. **Se placer dans le dossier `bin/`** :
4. Lancer l'exécutable : `projet_bioinfo` avec la commande `./projet_bioinfo`


### Paramétrage des variables globales

Le programme principal peut être configuré via des **variables globales** définies dans le fichier **`projet.h`** :

- `LONGUEUR_K_UPLET` : Définit la taille des k-uplets analysés.
- `LONGUEUR_MIN_MOTIF` : Longueur minimale d’un motif retenu.
- `X_HIT_AVANT_ENTRE` : Nombre minimal d'occurrences dans le génome complet (en excluant la séquence query en amont du gène) pour considérer un motif valide.
- `LONGUEUR_REGION_ETUIE_MOTIFS` : Longueur de la région promotrice analysée.
- `IDENTITE_MIN` : Pourcentage minimal d'identité pour détecter un gène.
- `BOITE_35` : Séquence de la boîte consensus -35.
- `BOITE_10` : Séquence de la boîte consensus -10.

Ces valeurs peuvent être modifiées en éditant directement **`projet.h`** avant la compilation.


### Dépendances

Pour compiler et exécuter le projet, vous aurez besoin des éléments suivants :

- Un compilateur C compatible, tel que GCC.
- Les librairies standard du langage C :
- `<stdio.h>`
- `<stdlib.h>`
- `<time.h>`
- `<unistd.h>`
- `<limits.h>`
- `<string.h>`
- `<stdbool.h>`
- Les fichiers suivants doivent être placés dans le dossier **`data/`** :
- Un fichier FASTA contenant le génome complet (`sequence_reelle.fna`).
- Un fichier FASTA contenant le gène d'intérêt (`gene.fna`).

