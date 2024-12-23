// recherche_consensus_box.c
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "projet.h"

/// Compare une sous-séquence à une boîte consensus avec une tolérance d'une erreur.
bool comparer_avec_marge_erreur(const char* sous_sequence, const char* sequence_consensus, int taille_consensus) {
    int erreurs = 0;
    for (int i = 0; i < taille_consensus; i++) {
        if (sous_sequence[i] != sequence_consensus[i]) {
            erreurs++;
        }
        if (erreurs > 1) {  // Plus d'une différence : la séquence ne correspond pas
            return false;
        }
    }
    return true;  // La séquence correspond avec au plus une erreur
}

// Fonction pour rechercher des boîtes consensus dans tout le génome ou en amont d'un gène
// Si base_amont == 0 la recherche est faite sur tout le génome
int rechercher_seq_consensus(const char* sequence_genome, const char* boite_35, const char* boite_10, int position_gene, int bases_amont) {
    int longueur_genome = strlen(sequence_genome);
    int longueur_boite_35 = strlen(boite_35);
    int longueur_boite_10 = strlen(boite_10);
    bool sequence_trouvee = false;
    int compteur = 0;  // Initialisation du compteur

    // Limiter la recherche à une région en amont si `bases_amont` > 0
    int debut_recherche = 0;
    int fin_recherche = longueur_genome;

    if (bases_amont > 0) {
        debut_recherche = (position_gene > bases_amont) ? position_gene - bases_amont : 0;
        fin_recherche = position_gene;
    }

    // Parcours de la séquence pour chercher les occurrences de la boîte -35
    for (int i = debut_recherche; i <= fin_recherche - longueur_boite_35; i++) {
        const char* sous_sequence_35 = sequence_genome + i;

        // Si une boîte -35 est trouvée (tolérance de 1 erreur)
        if (comparer_avec_marge_erreur(sous_sequence_35, boite_35, longueur_boite_35)) {

            // Recherche de la boîte -10 avec une distance de 15 à 19 nucléotides
            for (int distance = 15; distance <= 19; distance++) {
                int position_boite_10 = i + longueur_boite_35 + distance;

                // Vérifie que la boîte -10 ne dépasse pas la longueur du génome
                if (position_boite_10 + longueur_boite_10 > fin_recherche) {
                    continue;  // Ignore cette itération si la boîte -10 est hors limites
                }

                const char* sous_sequence_10 = sequence_genome + position_boite_10;

                // Si une boîte -10 est trouvée (tolérance de 1 erreur)
                if (comparer_avec_marge_erreur(sous_sequence_10, boite_10, longueur_boite_10)) {
                    sequence_trouvee = true;  // Une séquence valide a été trouvée
                    compteur++;  // Incrémenter le compteur pour chaque correspondance trouvée

                    // Ajuster les positions pour les rendre relatives au génome complet
                    int position_boite_35_absolue = i;
                    int position_boite_10_absolue = position_boite_10;

                    // Afficher les séquences et les positions trouvées
                    printf("Boîte -35 trouvée à la position %d : %.6s\n", position_boite_35_absolue, sous_sequence_35);
                    printf("Boîte -10 trouvée à la position %d : %.6s\n", position_boite_10_absolue, sous_sequence_10);
                }
            }
        }
    }

    // Si aucune séquence consensus n'a été trouvée, retourne -1
    if (!sequence_trouvee) {
        printf("Séquence consensus non trouvée\n");
        return -1;
    }

    // Si une ou plusieurs séquences ont été trouvées, retourne le nombre de séquences consensus trouvées
    printf("Nombre de séquences consensus trouvées : %d\n", compteur);
    return compteur;
}