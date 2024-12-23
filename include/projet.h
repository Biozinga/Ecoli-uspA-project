#ifndef PROJET_H
#define PROJET_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

// Constantes Générales
#define FOLD_CHANGE_MIN 1.5       ///< Seuil minimum du fold change pour l'enregistrement du motif dans le tableau des motifs retenu
#define LONGUEUR_K_UPLET 6      ///< Longueur des k-uplets lors du découpage de la séquence pour l'analyse des motifs
#define LONGUEUR_MIN_MOTIF 18    ///< Taille minimale devant être atteinte par un motif avant le test statistique

#define X_HIT_AVANT_ENTRE 20    ///< Nombre minimum d'occurrences pour qu'un hit soit retenu (même après extension)
#define LONGEUR_REGION_ETUIE_MOTIFS 1000 ///< Longueur de la région découpée en k-uplets (en amont du gène d'intérêt)
#define IDENTITE_MIN 0.90       ///< Identité minimale pour qu'un gène soit repéré dans la séquence

#define LONGEUR_SEQUENCE_ETUDIE_CONSENSUS 200  ///< Longueur de la région (en amont du gène d'intérêt) qui sera analysée pour trouver des séquences consensus
#define BOITE_35 "TTGACA"       ///< Séquence consensus (-35) de l'organisme étudié
#define BOITE_10 "TATAAT"       ///< Séquence consensus (-10) de l'organisme étudié

// Définition de la structure pour stocker un motif étendu
/**
 * @struct Motif
 * @brief Structure pour stocker les informations d'un motif étendu.
 */
typedef struct {
    char *sequence;                ///< Le motif étendu
    int start_pos;                 ///< Position du premier nucléotide du motif dans la séquence complète
    int length;                    ///< Longueur du motif
    double fold_change;            ///< Valeur du fold change
    int occurrences_reelles;       ///< Nombre d'occurrences réelles du motif
    int occurrences_aleatoires;    ///< Nombre d'occurrences dans la séquence aléatoire
} Motif;

// Prototypes des fonctions :

// Génération de séquence aléatoire

/**
 * @brief Génère un nucléotide aléatoire ('A', 'T', 'C' ou 'G').
 * @return Un caractère représentant un nucléotide aléatoire.
 */
char generer_nucleotide();

/**
 * @brief Génère une séquence aléatoire de nucléotides de longueur spécifiée.
 * @param longueur_seq La longueur de la séquence à générer.
 * @return Un pointeur vers la séquence de nucléotides générée (doit être libéré après utilisation).
 */
char* generer_sequence(int longueur_seq);

/**
 * @brief Écrit une séquence de nucléotides dans un fichier spécifié.
 * @param sequence La séquence de nucléotides à écrire.
 * @param chemin_fichier Le chemin vers le fichier de sortie.
 */
void ecrire_fichier_sequence(const char *sequence, const char *chemin_fichier);

// Analyse des motifs

/**
 * @brief Lit le contenu d'un fichier et renvoie une chaîne de caractères contenant les données.
 * @param chemin_fichier Le chemin vers le fichier à lire.
 * @return Un pointeur vers la chaîne contenant le contenu du fichier (doit être libéré après utilisation).
 */
char* lire_fichier(const char* chemin_fichier);

/**
 * @brief Recherche rapidement le nombre d'occurrences d'un motif dans une séquence complète, en excluant une sous-séquence spécifiée.
 * @param motif Le motif à rechercher.
 * @param sequence_complete La séquence complète dans laquelle rechercher.
 * @param start_sous_seq La position de début de la sous-séquence à exclure.
 * @param longueur_sous_sequence La longueur de la sous-séquence à exclure.
 * @return Le nombre d'occurrences du motif dans la séquence complète, hors de la sous-séquence exclue.
 */
int rechercher_motif_rapide(char *motif, char *sequence_complete, int start_sous_seq, int longueur_sous_sequence);

/**
 * @brief Calcule le fold change entre les occurrences réelles et aléatoires d'un motif.
 * @param occurrences_reelles Le nombre d'occurrences du motif dans la séquence réelle.
 * @param occurrences_aleatoires Le nombre d'occurrences du motif dans la séquence aléatoire.
 * @return La valeur du fold change.
 */
double calculer_fold_change(int occurrences_reelles, int occurrences_aleatoires);

/**
 * @brief Tente d'étendre un k-uplet à gauche et/ou à droite dans la séquence complète.
 * @param sequence_complete La séquence complète dans laquelle étendre le k-uplet.
 * @param longueur_total_genome La longueur totale du génome (séquence complète).
 * @param motif Un pointeur vers le motif à étendre.
 * @param occurrences_reelles Un pointeur vers le nombre d'occurrences réelles du motif (sera mis à jour).
 * @return 1 si une extension a été faite, 0 si aucune extension n'est possible, -1 en cas d'erreur.
 */
int etendre_k_uplet(char *sequence_complete, int longueur_total_genome, Motif *motif, int *occurrences_reelles);

/**
 * @brief Fonction de comparaison pour trier les motifs par position de départ.
 * @param a Pointeur vers le premier motif à comparer.
 * @param b Pointeur vers le second motif à comparer.
 * @return Une valeur négative si a < b, zéro si a == b, une valeur positive si a > b.
 */
int comparer_motifs(const void *a, const void *b);

/**
 * @brief Affiche une interface graphique résumant les motifs potentiels trouvés.
 * @param motifs_potentiels Un tableau de motifs potentiels.
 * @param nombre_motifs Le nombre de motifs dans le tableau.
 * @param position_gene La position du gène d'intérêt dans la séquence.
 */
void afficher_interface_graphique(Motif* motifs_potentiels, int nombre_motifs, int position_gene);

/**
 * @brief Traite les k-uplets dans une région de la séquence pour identifier des motifs potentiels.
 * @param sous_sequence La sous-séquence dans laquelle traiter les k-uplets.
 * @param longueur_sous_sequence La longueur de la sous-séquence.
 * @param sequence_complete La séquence complète (génome entier).
 * @param longueur_total_genome La longueur totale du génome.
 * @param longueur_k_uplet La longueur des k-uplets à utiliser.
 * @param longueur_min_motif La longueur minimale pour qu'un motif soit considéré.
 * @param sequence_aleatoire Une séquence aléatoire pour comparaison.
 * @param position_gene La position du gène d'intérêt dans la séquence complète.
 */
void traiter_k_uplets(char* sous_sequence, int longueur_sous_sequence, char* sequence_complete,
                      int longueur_total_genome, int longueur_k_uplet, int longueur_min_motif,
                      char* sequence_aleatoire, int position_gene);

/**
 * @brief Affiche une barre de progression dans la console.
 * @param current La valeur actuelle de progression.
 * @param total La valeur totale correspondant à 100% de progression.
 */
void afficher_progression(int current, int total);

// Recherche de gène

/**
 * @brief Calcule le pourcentage d'identité entre deux séquences de même longueur.
 * @param seq1 La première séquence.
 * @param seq2 La seconde séquence.
 * @param longueur La longueur des séquences à comparer.
 * @return Le pourcentage d'identité entre les deux séquences.
 */
float calculer_identite(const char* seq1, const char *seq2, int longueur);

/**
 * @brief Recherche la position d'un gène dans une séquence génomique en fonction d'une identité minimale.
 * @param sequence_genome La séquence génomique complète.
 * @param sequence_gene La séquence du gène à rechercher.
 * @param identite_min Le pourcentage d'identité minimale requis pour considérer une correspondance.
 * @return La position du gène dans la séquence génomique, ou -1 si non trouvé ou trouvé plusieurs fois.
 */
int rechercher_gene(const char* sequence_genome, const char *sequence_gene, double identite_min);

// Recherche de séquence consensus

/**
 * @brief Compare une sous-séquence avec une séquence consensus avec une marge d'erreur de 1 base.
 * @param sous_sequence La sous-séquence extraite du génome.
 * @param sequence_consensus La séquence consensus à comparer.
 * @param taille_consensus La longueur de la séquence consensus.
 * @return true si les séquences correspondent avec au plus 1 différence, false sinon.
 */
bool comparer_avec_marge_erreur(const char* sous_sequence, const char* sequence_consensus, int taille_consensus);

/**
 * @brief Recherche les séquences consensus de type boîte -35 et boîte -10 en amont d'un gène.
 * @param sequence_genome La séquence génomique complète.
 * @param boite_35 La séquence consensus de la boîte -35.
 * @param boite_10 La séquence consensus de la boîte -10.
 * @param position_gene La position du gène dans la séquence génomique.
 * @param bases_amont Le nombre de bases en amont du gène à analyser.
 * @return Le nombre de paires de boîtes consensus trouvées, ou -1 si aucune n'est trouvée.
 */
int rechercher_seq_consensus(const char* sequence_genome, const char* boite_35, const char* boite_10,
                             int position_gene, int bases_amont);

// Tests des fonctions :

// Tests de séquences consensus

/**
 * @brief Teste la détection de boîtes consensus avec une distance de 15 nucléotides entre les boîtes -35 et -10.
 */
void test_seq_consensus_distance_15();

/**
 * @brief Teste la détection de boîtes consensus avec une distance de 19 nucléotides entre les boîtes -35 et -10.
 */
void test_seq_consensus_distance_19();

/**
 * @brief Teste la détection de boîtes consensus chevauchantes avec des distances variables.
 */
void test_seq_consensus_chevauchante();

/**
 * @brief Teste la robustesse de la détection de boîtes consensus dans des séquences riches en boîtes -35 et -10.
 */
void test_seq_consensus_critique_seq();

/**
 * @brief Teste la recherche de boîtes consensus en amont avec une limitation du nombre de bases à analyser.
 */
void test_rechercher_seq_consensus_en_amont();

/**
 * @brief Teste la recherche de boîtes consensus avec un gène chevauchant la fin de la boîte -10.
 */
void test_rechercher_seq_consensus_en_amont_critique();

// Tests de recherche de gène

/**
 * @brief Teste la détection d'un gène présent avec une identité parfaite.
 * @return 0 si le test est réussi.
 */
int test_presence_gene();

/**
 * @brief Teste la détection d'un gène présent plusieurs fois dans le génome.
 * @return 0 si le test est réussi.
 */
int test_presence_multiple_gene();

/**
 * @brief Teste l'absence de détection d'un gène lorsque l'identité est insuffisante.
 * @return 0 si le test est réussi.
 */
int test_abscence_gene();

/**
 * @brief Teste l'absence de détection d'un gène présent mais avec une identité insuffisante (cas critique).
 * @return 0 si le test est réussi.
 */
int test_abscence_gene_critique();

// Tests d'analyse de motifs

/**
 * @brief Teste la fonction de recherche rapide de motifs dans une séquence.
 */
void test_rechercher_motif_rapide();

/**
 * @brief Teste la fonction de calcul du fold change entre les occurrences réelles et aléatoires.
 */
void test_calculer_fold_change();

/**
 * @brief Teste la fonction d'extension et de traitement d'un k-uplet.
 */
void test_etendre_et_traiter_k_uplets();

#endif // PROJET_H
