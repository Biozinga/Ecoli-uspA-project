// recherche_motifs.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>

#include "projet.h"

#define LONGEUR_REGION_ETUIE_MOTIFS 1000 // Longueur de la région découpée en k-uplets (en amont du gène d'intérêt)

// Fonction pour lire un fichier FASTA ou texte et en extraire la séquence d'ADN
char* lire_fichier(const char* chemin_fichier) {
    FILE *fichier = fopen(chemin_fichier, "r");
    if (fichier == NULL) {
        perror("Erreur lors de l'ouverture du fichier");
        return NULL;
    }

    fseek(fichier, 0, SEEK_END);
    long taille_fichier = ftell(fichier);
    rewind(fichier);

    char *sequence = (char*)malloc((taille_fichier + 1) * sizeof(char));
    if (sequence == NULL) {
        printf("Erreur d'allocation mémoire\n");
        fclose(fichier);
        return NULL;
    }

    char ligne[256];
    int index = 0;
    while (fgets(ligne, sizeof(ligne), fichier)) {
        if (ligne[0] == '>') continue;  // Ignorer les lignes de description (FASTA)
        for (int i = 0; ligne[i] != '\0'; i++) {
            if (ligne[i] != '\n' && ligne[i] != ' ') {
                char nucleotide = toupper(ligne[i]);
                if (nucleotide == 'A' || nucleotide == 'C' || nucleotide == 'G' || nucleotide == 'T') {
                    sequence[index++] = nucleotide;
                }
            }
        }
    }

    sequence[index] = '\0';
    fclose(fichier);
    return sequence;
}

// Fonction modifiée pour rechercher un motif dans une séquence en excluant une partie du génome (k-uplets)
int rechercher_motif_rapide(char *motif, char *sequence_complete, int start_sous_seq, int longueur_sous_sequence) {
    int occurrences = 0;
    char *position = sequence_complete;

    // Boucle de recherche pour le motif dans la séquence complète
    while ((position = strstr(position, motif)) != NULL) {
        int position_index = position - sequence_complete;  // Calculer la position de l'occurrence trouvée
        
        // Vérifier si l'occurrence est dans la région à exclure (le k-uplet)
        if (position_index >= start_sous_seq && position_index < start_sous_seq + longueur_sous_sequence) {
            // Si l'occurrence est dans la sous-séquence découpée (k-uplets), on l'ignore
            position++;  // Continuer la recherche après cette occurrence
            continue;
        }

        // Si l'occurrence n'est pas dans la sous-séquence à exclure, on la compte
        occurrences++;
        position++;  // Continuer la recherche après cette occurrence
    }

    return occurrences;
}

// Fonction pour calculer le fold change
double calculer_fold_change(int occurrences_reelles, int occurrences_aleatoires) {
    double epsilon = 1e-6;
    return (double)occurrences_reelles / (occurrences_aleatoires + epsilon); // epsilon pour evité la division par 0
}

// Fonction pour essayer d'étendre un k-uplet à gauche et à droite
int etendre_k_uplet(char *sequence_complete, int longueur_total_genome, Motif *motif, int *occurrences_reelles) {
    int left = motif->start_pos;         // Position pour étendre à gauche
    int right = motif->start_pos + motif->length;  // Position pour étendre à droite
    int extended = 0;  // Flag pour savoir si une extension a été faite
    char *sequence_temp = NULL;
    int new_occurrences = 0;

    // Tenter d'étendre à gauche
    if (left > 0) {
        int new_length = motif->length + 1;  // Étendre d'une base à gauche
        sequence_temp = malloc((new_length + 1) * sizeof(char));
        if (sequence_temp == NULL) {
            printf("Erreur : échec de l'allocation mémoire\n");
            return -1;
        }

        // Ajouter la base à gauche
        sequence_temp[0] = sequence_complete[left - 1];
        memcpy(sequence_temp + 1, motif->sequence, motif->length);
        sequence_temp[new_length] = '\0';

        // Rechercher dans la séquence complète
        new_occurrences = rechercher_motif_rapide(sequence_temp, sequence_complete, motif->start_pos, motif->length);

        if (new_occurrences >= X_HIT_AVANT_ENTRE) {
            free(motif->sequence);
            motif->sequence = sequence_temp;
            motif->start_pos--;  // Mettre à jour la position de départ
            motif->length = new_length;
            *occurrences_reelles = new_occurrences;
            extended = 1;  // Extension réussie
        } else {
            free(sequence_temp);  // Libérer la mémoire si échec
        }
    }

    // Tenter d'étendre à droite
    if (right < longueur_total_genome && !extended) {  // Étendre à droite seulement si pas étendu à gauche
        int new_length = motif->length + 1;  // Étendre d'une base à droite
        sequence_temp = malloc((new_length + 1) * sizeof(char));
        if (sequence_temp == NULL) {
            printf("Erreur : échec de l'allocation mémoire\n");
            return -1;
        }

        // Copier l'ancien motif et ajouter la base à droite
        memcpy(sequence_temp, motif->sequence, motif->length);
        sequence_temp[new_length - 1] = sequence_complete[right];
        sequence_temp[new_length] = '\0';

        // Rechercher dans la séquence complète
        new_occurrences = rechercher_motif_rapide(sequence_temp, sequence_complete, motif->start_pos, motif->length);

        if (new_occurrences >= X_HIT_AVANT_ENTRE) {
            free(motif->sequence);
            motif->sequence = sequence_temp;
            motif->length = new_length;
            *occurrences_reelles = new_occurrences;
            extended = 1;  // Extension réussie
        } else {
            free(sequence_temp);  // Libérer la mémoire si échec
        }
    }

    return extended;
}


// Fonction de comparaison pour le tri des motifs par position
int comparer_motifs(const void *a, const void *b) {
    Motif *motifA = (Motif *)a;
    Motif *motifB = (Motif *)b;
    return motifA->start_pos - motifB->start_pos;
}

// Fonction pour afficher la barre de progression
void afficher_progression(int current, int total) {
    int largeur_barre = 50; // Largeur de la barre de progression
    float progression = (float)current / total;
    int position = (int)(progression * largeur_barre);

    printf("\r["); // Retour au début de la ligne
    for (int i = 0; i < largeur_barre; i++) {
        if (i < position) {
            printf("=");  // Simple '=' pour la barre
        } else {
            printf(" ");  // Espaces pour combler la barre
        }
    }
    printf("] %d%%", (int)(progression * 100));
    fflush(stdout);
}

// Fonction pour afficher l'interface graphique des résultats
void afficher_interface_graphique(Motif* motifs_potentiels, int nombre_motifs, int position_gene) {
    if (nombre_motifs == 0) {
        printf("Aucun motif potentiel trouvé.\n");
        return;
    }

    // Trouver le motif avec le fold change le plus élevé
    Motif* meilleur_fold_change = &motifs_potentiels[0];
    Motif* motif_plus_proche_gene = &motifs_potentiels[0];
    Motif* motif_le_plus_long = &motifs_potentiels[0];

    for (int i = 1; i < nombre_motifs; i++) {
        if (motifs_potentiels[i].fold_change > meilleur_fold_change->fold_change) {
            meilleur_fold_change = &motifs_potentiels[i];
        }
        // Calculer la distance au gène
        int distance_courante = abs(motifs_potentiels[i].start_pos - position_gene);
        int distance_minimale = abs(motif_plus_proche_gene->start_pos - position_gene);
        if (distance_courante < distance_minimale) {
            motif_plus_proche_gene = &motifs_potentiels[i];
        }
        if (motifs_potentiels[i].length > motif_le_plus_long->length) {
            motif_le_plus_long = &motifs_potentiels[i];
        }
    }

    // Affichage graphique final avec doubles lignes '======'
    printf("\n====== Motif avec le taux de fold change le plus élevé ======\n");
    printf("Motif : %s\n", meilleur_fold_change->sequence);
    printf("Position : %d\n", meilleur_fold_change->start_pos);
    printf("Longueur : %d\n", meilleur_fold_change->length);
    printf("Fold Change : %.2f\n", meilleur_fold_change->fold_change);
    printf("\n");

    printf("====== Motif le plus proche du gène ======\n");
    printf("Motif : %s\n", motif_plus_proche_gene->sequence);
    printf("Position : %d\n", motif_plus_proche_gene->start_pos);
    printf("Longueur : %d\n", motif_plus_proche_gene->length);
    printf("Fold Change : %.2f\n", motif_plus_proche_gene->fold_change);
    printf("\n");

    printf("====== Motif ayant la longueur la plus élevée ======\n");
    printf("Motif : %s\n", motif_le_plus_long->sequence);
    printf("Position : %d\n", motif_le_plus_long->start_pos);
    printf("Longueur : %d\n", motif_le_plus_long->length);
    printf("Fold Change : %.2f\n", motif_le_plus_long->fold_change);
    printf("\n");
}

// Fonction pour traiter les k-uplets dans une région de la séquence
void traiter_k_uplets(char* sous_sequence, int longueur_sous_sequence, char* sequence_complete,
                      int longueur_total_genome, int longueur_k_uplet, int longueur_min_motif,
                      char* sequence_aleatoire, int position_gene) {
    // Calculer le début de la sous-séquence dans la séquence complète
    int start_sous_seq = (position_gene >= LONGEUR_REGION_ETUIE_MOTIFS) ? (position_gene - LONGEUR_REGION_ETUIE_MOTIFS) : 0;

    // Tableau dynamique pour stocker les motifs potentiels
    Motif *motifs_potentiels = NULL;
    int nombre_motifs = 0;
    int capacite_motifs = 0;

    int total_k_uplets = longueur_sous_sequence - longueur_k_uplet + 1;

    // Afficher un séparateur avant le début de l'analyse
    printf("Début de l'analyse des k-uplets\n");

    // Traiter les k-uplets dans la région de la séquence
    for (int i = 0; i <= longueur_sous_sequence - longueur_k_uplet; i++) {
        // Afficher la progression
        afficher_progression(i + 1, total_k_uplets);

        // Extraction du k-uplet
        char k_uplet[longueur_k_uplet + 1];
        strncpy(k_uplet, &sous_sequence[i], longueur_k_uplet);
        k_uplet[longueur_k_uplet] = '\0';

        // Définir les paramètres pour exclure cette région lors de la recherche
        int start_pos_exclusion = start_sous_seq + i;  // Position exacte du k-uplet dans la séquence complète
        int longueur_exclusion = longueur_k_uplet;

        // Rechercher ce k-uplet dans la séquence complète
        int occurrences_reelles = rechercher_motif_rapide(k_uplet, sequence_complete, start_pos_exclusion, longueur_exclusion);

        // Rechercher dans la séquence aléatoire
        int occurrences_aleatoires = rechercher_motif_rapide(k_uplet, sequence_aleatoire, 0, 0);

        // Si le motif est trouvé plusieurs fois dans la séquence réelle, étendre et tester
        if (occurrences_reelles >= X_HIT_AVANT_ENTRE) {
            // Calculer la position de départ dans le génome complet
            int start_pos_in_genome = start_sous_seq + i;

            // Créer un motif initial à partir du k-uplet
            Motif motif = {strdup(k_uplet), start_pos_in_genome, longueur_k_uplet, 0.0, occurrences_reelles, occurrences_aleatoires};

            int extension_possible = 1;

            // Essayer d'étendre le k-uplet à gauche et à droite jusqu'à ce qu'on ne puisse plus l'étendre
            while (extension_possible) {
                extension_possible = etendre_k_uplet(sequence_complete, longueur_total_genome, &motif, &occurrences_reelles);
                if (extension_possible == -1) {
                    // Erreur lors de l'extension
                    break;
                }
            }

            // Vérifier si le motif atteint la longueur minimale avant de faire le test statistique
            if (motif.length >= longueur_min_motif) {
                // Recalculer les occurrences réelles pour le motif étendu
                occurrences_reelles = rechercher_motif_rapide(motif.sequence, sequence_complete,  motif.start_pos, motif.length);
                motif.occurrences_reelles = occurrences_reelles;

                occurrences_aleatoires = rechercher_motif_rapide(motif.sequence, sequence_aleatoire, 0, 0);
                motif.occurrences_aleatoires = occurrences_aleatoires;

                double fold_change = calculer_fold_change(occurrences_reelles, occurrences_aleatoires);

                if (fold_change > FOLD_CHANGE_MIN) {
                    // Stocker le motif dans le tableau
                    if (nombre_motifs == capacite_motifs) {
                        // Augmenter la capacité du tableau
                        capacite_motifs = capacite_motifs == 0 ? 10 : capacite_motifs * 2;
                        motifs_potentiels = realloc(motifs_potentiels, capacite_motifs * sizeof(Motif));
                        if (motifs_potentiels == NULL) {
                            printf("Erreur d'allocation mémoire pour le tableau de motifs\n");
                            free(motif.sequence);
                            break;
                        }
                    }

                    // Ajouter le motif au tableau
                    motif.fold_change = fold_change;
                    motifs_potentiels[nombre_motifs++] = motif;
                } else {
                    // Libérer la mémoire si le motif n'est pas retenu
                    free(motif.sequence);
                }
            } else {
                free(motif.sequence);
            }
        }
    }

    // Ajouter une nouvelle ligne après la barre de progression
    printf("\n");

    // Afficher le résumé des motifs après tri par position
    if (nombre_motifs > 0) {
        printf("\n====== Analyse terminée ======\n");
        // Trier les motifs par position
        qsort(motifs_potentiels, nombre_motifs, sizeof(Motif), comparer_motifs);

        printf("Résumé des motifs potentiels retenus :\n");

        // Enregistrement des 3 meilleurs motifs dans motif_retenu.txt au forma FASTA
        FILE* fichier_motifs = fopen("../data/motif_retenu.txt", "w");  // AJOUT : Ouverture du fichier
        if (fichier_motifs == NULL) {
            perror("Impossible de créer ou d'ouvrir le fichier data/motif_retenu.txt");
        }

        // Affichage CONSOLE de tous les motifs 
        for (int i = 0; i < nombre_motifs; i++) {
            printf("Motif : %s\n", motifs_potentiels[i].sequence);
            printf("Position : %d\n", motifs_potentiels[i].start_pos);
            printf("Longueur : %d\n", motifs_potentiels[i].length);
            printf("Fold Change : %.2f\n", motifs_potentiels[i].fold_change);
            printf("Occurrences réelles : %d\n", motifs_potentiels[i].occurrences_reelles);
            printf("Occurrences aléatoires : %d\n", motifs_potentiels[i].occurrences_aleatoires);
            printf("-----------------------------\n");
        }

        // Recherche des 3 motifs « meilleurs » (mêmes calculs que dans afficher_interface_graphique)
        Motif* meilleur_fold_change = &motifs_potentiels[0];
        Motif* motif_plus_proche_gene = &motifs_potentiels[0];
        Motif* motif_le_plus_long = &motifs_potentiels[0];

        for (int i = 1; i < nombre_motifs; i++) {
            if (motifs_potentiels[i].fold_change > meilleur_fold_change->fold_change) {
                meilleur_fold_change = &motifs_potentiels[i];
            }
            int distance_courante = abs(motifs_potentiels[i].start_pos - position_gene);
            int distance_minimale = abs(motif_plus_proche_gene->start_pos - position_gene);
            if (distance_courante < distance_minimale) {
                motif_plus_proche_gene = &motifs_potentiels[i];
            }
            if (motifs_potentiels[i].length > motif_le_plus_long->length) {
                motif_le_plus_long = &motifs_potentiels[i];
            }
        }


        if (fichier_motifs != NULL) {
            // Motif avec le plus fort fold change
            fprintf(fichier_motifs, ">Motif avec le taux de fold change le plus élevé\n");
            fprintf(fichier_motifs, "%s\n\n", meilleur_fold_change->sequence);

            // Motif le plus proche du gène
            fprintf(fichier_motifs, ">Motif le plus proche du gène\n");
            fprintf(fichier_motifs, "%s\n\n", motif_plus_proche_gene->sequence);

            // Motif ayant la longueur la plus élevée
            fprintf(fichier_motifs, ">Motif ayant la longueur la plus élevée\n");
            fprintf(fichier_motifs, "%s\n\n", motif_le_plus_long->sequence);

            fclose(fichier_motifs);  // AJOUT : Fermeture du fichier
        }


        // Afficher l'interface graphique finale
        afficher_interface_graphique(motifs_potentiels, nombre_motifs, position_gene);

        // Libérer la mémoire des motifs
        for (int i = 0; i < nombre_motifs; i++) {
            free(motifs_potentiels[i].sequence);
        }
        free(motifs_potentiels);
    } else {
        printf("\n====== Analyse terminée ======\n");
        printf("Aucun motif potentiel retenu.\n");
    }
}
