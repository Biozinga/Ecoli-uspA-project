#include "projet.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <limits.h>
#include <string.h>

int main() {
    // ==================== TEST DES FONCTIONS ====================
    printf("\n========================================================================================\n");
    printf("         Test des fonctions avant lancement du programme\n");
    printf("========================================================================================\n\n");

    // Tests des fonctions de recherche de séquences consensus
    test_seq_consensus_distance_15();
    test_seq_consensus_distance_19();
    test_seq_consensus_chevauchante();
    test_seq_consensus_critique_seq();
    test_rechercher_seq_consensus_en_amont();
    test_rechercher_seq_consensus_en_amont_critique();

    // Tests des fonctions de recherche de gènes
    test_presence_gene();
    test_presence_multiple_gene();
    test_abscence_gene();
    test_abscence_gene_critique();

    // Tests des fonctions de recherche de motifs
    test_rechercher_motif_rapide();
    test_calculer_fold_change();
    
    // ATTENTION, ces tests arriveront dans la prochaine version du programme.
    
    // Tests des fonctions d'extention  et de traitement des k-uplets
    //test_etendre_et_traiter_k_uplets();

    // Message final indiquant que tous les tests sont réussis
    printf("\n========================================================================================\n");
    printf("   Tous les tests ont été passés avec succès\n");
    printf("========================================================================================\n\n");

    // ==================== PROGRAMME PRINCIPAL ====================
    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) == NULL) {
        perror("Erreur lors de l'obtention du répertoire de travail");
        return 1;
    }

    // Chemins des fichiers
    char chemin_sequence_reelle[PATH_MAX];
    char chemin_sequence_aleatoire[PATH_MAX];
    char chemin_gene[PATH_MAX];

    snprintf(chemin_sequence_reelle, PATH_MAX, "%s/../data/sequence_reelle.fna", cwd);
    snprintf(chemin_sequence_aleatoire, PATH_MAX, "%s/../data/sequence_aleatoire.fna", cwd);
    snprintf(chemin_gene, PATH_MAX, "%s/../data/gene.fna", cwd);

    // Initialisation de la graine aléatoire
    srand(time(NULL));

    // Génération de la séquence aléatoire
    // Lecture du génome complet étudié
    char* sequence_reelle = lire_fichier(chemin_sequence_reelle);
    if (!sequence_reelle) {
        printf("Erreur lors de la lecture de la séquence réelle.\n");
        return 1;
    }
    int longueur_sequence_aleatoire = strlen(sequence_reelle);
    char *sequence = generer_sequence(longueur_sequence_aleatoire);
    if (sequence == NULL) {
        free(sequence_reelle);
        return 1;
    }

    ecrire_fichier_sequence(sequence, chemin_sequence_aleatoire);
    free(sequence);

    // Lecture des séquences aléatoire et du gène d'intérêt
    char* sequence_alea = lire_fichier(chemin_sequence_aleatoire);
    char* sequence_gene = lire_fichier(chemin_gene);

    if (!sequence_alea || !sequence_gene) {
        printf("Erreur lors de la lecture des fichiers de séquence.\n");
        free(sequence_reelle);
        free(sequence_alea);
        free(sequence_gene);
        return 1;
    }

    // Séparateur
    printf("\n========================================================================================\n");
    printf("   Recherche du gène dans la séquence réelle d'E. coli\n");
    printf("========================================================================================\n\n");

    // Paramètres pour la recherche du gène
    double identite_minimale = IDENTITE_MIN;  // Pourcentage d'identité minimale

    // Recherche du gène dans la séquence réelle et enregistrement de sa position
    int position_gene = rechercher_gene(sequence_reelle, sequence_gene, identite_minimale);

    // Vérifier que la position du gène est valide
    int longueur_sequence_reelle = strlen(sequence_reelle);
    if (position_gene < 0 || position_gene > longueur_sequence_reelle) {
        printf("Erreur : position du gène invalide.\n");
        free(sequence_reelle);
        free(sequence_alea);
        free(sequence_gene);
        return 1;
    }

    // === Recherche de séquences consensus avant l'analyse des motifs ===
    printf("\n========================================================================================\n");
    printf("   Recherche de la présence de boîtes consensus en amont du gène\n");
    printf("========================================================================================\n\n");

    // Paramètres pour la recherche de séquences consensus
    int bases_amont_consensus = LONGEUR_SEQUENCE_ETUDIE_CONSENSUS;  // Nombre de bases en amont pour la recherche
    char* boite_35 = BOITE_35;  // Séquence de la boîte -35
    char* boite_10 = BOITE_10;  // Séquence de la boîte -10

    // Appel à la fonction rechercher_seq_consensus
    int nombre_sequences = rechercher_seq_consensus(sequence_reelle, boite_35, boite_10, position_gene, bases_amont_consensus);

    if (nombre_sequences == -1) {
        printf("Aucune séquence consensus trouvée en amont du gène.\n");
    } else {
        printf("Nombre total de séquences consensus trouvées : %d\n", nombre_sequences);
    }

    // Séparateur pour l'analyse des motifs
    printf("\n========================================================================================\n");
    printf("   Analyse des motifs dans la région promotrice\n");
    printf("========================================================================================\n\n");

    // Vérifier que la région en amont ne dépasse pas les limites de la séquence
    int taille_region = LONGEUR_REGION_ETUIE_MOTIFS;  // Taille de la région en amont à analyser
    if (position_gene - taille_region < 0) {
        printf("Erreur : la région en amont dépasse les limites de la séquence\n");
        free(sequence_reelle);
        free(sequence_alea);
        free(sequence_gene);
        return 1;
    }

    // Extraire la sous-séquence à analyser en amont du gène
    char* sous_sequence = (char*)malloc((taille_region + 1) * sizeof(char));
    if (sous_sequence == NULL) {
        printf("Erreur d'allocation mémoire pour la sous-séquence\n");
        free(sequence_reelle);
        free(sequence_alea);
        free(sequence_gene);
        return 1;
    }

    strncpy(sous_sequence, &sequence_reelle[position_gene - taille_region], taille_region);
    sous_sequence[taille_region] = '\0';

    // Paramètres pour l'analyse des k-uplets
    int longueur_k_uplet = LONGUEUR_K_UPLET;     // Longueur des k-uplets à rechercher
    int longueur_min_motif = LONGUEUR_MIN_MOTIF; // Longueur minimale d'un motif pour le test statistique

    // Appeler la fonction pour traiter les k-uplets en excluant la région spécifique
    traiter_k_uplets(sous_sequence, taille_region, sequence_reelle,
                     longueur_sequence_reelle, longueur_k_uplet, longueur_min_motif,
                     sequence_alea, position_gene);

    // Libérer la mémoire
    free(sous_sequence);

    // === Fin de l'intégration du module analyse_motifs ===

    printf("\nProgramme terminé avec succès.\n");

    // Libération de la mémoire
    free(sequence_reelle);
    free(sequence_alea);
    free(sequence_gene);

    return 0;
}
