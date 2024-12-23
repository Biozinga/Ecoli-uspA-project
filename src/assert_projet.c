// assert_projet.c
#include <stdio.h>
#include <assert.h>
#include "projet.h"  

// assert recherche_consensus_box.c

// Test d'une séquence avec une distance de 15 nucléotides entre les boîtes -35 et -10
void test_seq_consensus_distance_15() {
    const char* genome = "TTGACATACGGGAGGGGGGGGTATAATGCGTAC";  // Distance de 15
    const char* boite_35 = "TTGACA";  // Séquence boîte -35
    const char* boite_10 = "TATAAT";  // Séquence boîte -10
    int position_gene = strlen(genome) - 1;  // Position du dernier caractère du gène

    // Appel de la fonction avec `bases_amont` à 0 (aucune limitation de la recherche)
    assert(rechercher_seq_consensus(genome, boite_35, boite_10, position_gene, 0) == 1);
}

// Test d'une séquence avec une distance de 19 nucléotides entre les boîtes -35 et -10
void test_seq_consensus_distance_19() {
    const char* genome = "TTGACATCGATCGAGAGTAAGATCATATAATGCGTAC";  // Distance de 19
    const char* boite_35 = "TTGACA";  // Séquence boîte -35
    const char* boite_10 = "TATAAT";  // Séquence boîte -10
    int position_gene = strlen(genome) - 1;  // Position du dernier caractère du gène

    // Appel de la fonction avec `bases_amont` à 0
    assert(rechercher_seq_consensus(genome, boite_35, boite_10, position_gene, 0) == 1);
}

// Test de séquences chevauchantes (avec un cadre de lecture décalé) avec distance de 15 et 19 nucléotides
void test_seq_consensus_chevauchante() {
    const char* genome = "TTGACATATTGACATTTGACACGGTATAATATAATTGC";  // Séquences chevauchantes
    const char* boite_35 = "TTGACA";  // Séquence boîte -35
    const char* boite_10 = "TATAAT";  // Séquence boîte -10
    int position_gene = strlen(genome) - 1;  // Position du dernier caractère du gène

    // Appel de la fonction avec `bases_amont` à 0
    assert(rechercher_seq_consensus(genome, boite_35, boite_10, position_gene, 0) == 2);  
}

// Test critique pour vérifier la robustesse du code avec des séquences riches en boîtes -35 et -10, mais dont la distance ne permet pas la détection.
void test_seq_consensus_critique_seq() {
    const char* genome = "TTGACATATAATTTGACATATAATTTGACATATAATTTGACATATAATGCGGATATAAT";  // Génome très riche en -35 -10 
    const char* boite_35 = "TTGACA";  // Séquence boîte -35
    const char* boite_10 = "TATAAT";  // Séquence boîte -10
    int position_gene = strlen(genome) - 1;  // Position du dernier caractère du gène

    // Appel de la fonction avec `bases_amont` à 0
    assert(rechercher_seq_consensus(genome, boite_35, boite_10, position_gene, 0) == -1);  
}

// Fonction de test pour la recherche de boîtes consensus avec de petites séquences
void test_rechercher_seq_consensus_en_amont() {
    const char* genome = "AATTGACACACCGGCATTACTTAAAGTATAATGCC";  // Génome simple avec une boîte -35 et une boîte -10 mais base en amont 
    const char* boite_35 = "TTGACA";  // Boîte -35                  trop restrictive
    const char* boite_10 = "TATAAT";  // Boîte -10
    int position_gene = 32;  // Position du gène
    int bases_amont = 20;  // Limite de 20 bases en amont

    // Appel de la fonction avec une limite de 20 bases en amont
    assert(rechercher_seq_consensus(genome, boite_35, boite_10, position_gene, bases_amont) == -1);
}

// Fonction de test pour la recherche de boîtes consensus avec gene chevauchant la fin de boite -10
void test_rechercher_seq_consensus_en_amont_critique() {
    const char* genome = "AATTGACACACCGGCATTACTTAAAGTATAATGCC";  // Génome simple avec une boîte -35 et une boîte -10 mais base en amont 
    const char* boite_35 = "TTGACA";  // Boîte -35                  trop restrictive
    const char* boite_10 = "TATAAT";  // Boîte -10
    int position_gene = 31;  // Position du gène sur fin de boite
    int bases_amont = 32;  // Recherche sur toute la sequence -1

    // Appel de la fonction avec une limite de 20 bases en amont
    assert(rechercher_seq_consensus(genome, boite_35, boite_10, position_gene, bases_amont) == -1);
}

// assert recherche_gene.c

// Utilisation des <= dans un contexte de float
// Test d'un gène parfaitement trouvé avec une identité minimale de 100%
int test_presence_gene() {
    const char* genome = "ACGTACGTCCCACTACGTACTACGTTAACGTCCAATCTATATGGATCATGTTACCCCATCGGAGTACTACGTATTTCGTACGT";
    const char* gene = "ACGTCCAATCTATATGGATCATGTTACCCCATCGGAGTA";
    assert(calculer_identite(genome + 27, gene, strlen(gene)) >= 0.999 && calculer_identite(genome + 28, gene, strlen(gene)) <= 1.001);
    assert(rechercher_gene(genome, gene, 1.0) == 28);  // Recherchera une correspondance parfaite
    // Gène présent en position 28 avec une identité de 100%
    return 0;
    
}

// Test d'un gène trouvé plusieurs fois
int test_presence_multiple_gene() {
    const char* genome = "ACGTACGTCCCACTAACGTCCAATCTATATGGATCATGTTACCCCATCGGAGTACGTACTACGTTAACGTCCAATCTATATGGATCATGTTACCCCATCGGAGTACTACGTATTTCGTACGT";
    const char* gene = "ACGTCCAATCTATATGGATCATGTTACCCCATCGGAGTA";
    assert(rechercher_gene(genome, gene, 1.0) == -1);  
    // Gène présent deux fois la fonction dois return -1
    return 0;
    
}

// Test d'un gène non trouvé car identité insuffisante
int test_abscence_gene() {
    const char* genome = "ACGTACGTACGT";
    const char* gene = "ACGGG";  // Une seule différence
    assert(calculer_identite(genome, gene, strlen(gene)) >= 0.599 && calculer_identite(genome, gene, strlen(gene)) <= 0.601);
    assert(rechercher_gene(genome, gene, 0.8) == -1);  
    // Pas d'assertion directe ici, le résultat doit être affiché par rechercher_gene
    // Le gène testé à 60% d'identité, donc ne doit pas être trouvé
    return 0;
}

// Test avec un gène qui est présent mais avec une identité insufisante (78.57%)
int test_abscence_gene_critique() {
    const char* genome = "ACGTACGTACGTTG";
    const char* gene = "ACGTATGAACGATG";  // 11 bases identiques sur 14 => 78.57% d'identité
    assert(calculer_identite(genome, gene, strlen(gene)) >= 0.784 && calculer_identite(genome, gene, strlen(gene)) <= 0.786);
    assert(rechercher_gene(genome, gene, 0.8) == -1);  // Ne doit pas être trouvé
    // Pas d'assertion directe ici, le résultat doit être affiché par rechercher_gene
    return 0;
}

// assert motif_analyse.c

// Test de la fonction rechercher_motif_rapide avec des printf pour le débogage
void test_rechercher_motif_rapide() {
    printf("=== Début du test de rechercher_motif_rapide ===\n");
    
    char *sequence = "ACGTACGTACGTACGT";
    char *motif = "ACGT";
    
    printf("Test du motif '%s' dans la séquence '%s'\n", motif, sequence);
    int occurrences = rechercher_motif_rapide(motif, sequence, 0, 0);
    printf("Occurrences trouvées : %d\n", occurrences);
    assert(occurrences == 4);  // Le motif "ACGT" apparaît 4 fois
    
    motif = "CGT";
    printf("Test du motif '%s' dans la séquence '%s'\n", motif, sequence);
    occurrences = rechercher_motif_rapide(motif, sequence, 0, 0);
    printf("Occurrences trouvées : %d\n", occurrences);
    assert(occurrences == 4);  // Le motif "CGT" apparaît 4 fois
    
    motif = "TACG";
    printf("Test du motif '%s' dans la séquence '%s'\n", motif, sequence);
    occurrences = rechercher_motif_rapide(motif, sequence, 0, 0);
    printf("Occurrences trouvées : %d\n", occurrences);
    assert(occurrences == 3);  // Le motif "TACG" apparaît 3 fois
    
    motif = "AAA";
    printf("Test du motif '%s' dans la séquence '%s'\n", motif, sequence);
    occurrences = rechercher_motif_rapide(motif, sequence, 0, 0);
    printf("Occurrences trouvées : %d\n", occurrences);
    assert(occurrences == 0);  // Le motif "AAA" n'apparaît pas
    
    printf("Test de rechercher_motif_rapide passé avec succès.\n");
}

// Test de la fonction calculer_fold_change avec des printf pour le débogage
void test_calculer_fold_change() {
    printf("=== Début du test de calculer_fold_change ===\n");
    
    int occurrences_reelles = 20;
    int occurrences_aleatoires = 10;
    double fold_change = calculer_fold_change(occurrences_reelles, occurrences_aleatoires);
    printf("Fold change pour %d occurrences réelles et %d occurrences aléatoires : %.2f\n", occurrences_reelles, occurrences_aleatoires, fold_change);
    assert(fold_change > 1.99 && fold_change < 2.01);  // Fold change ≈ 2.0
    
    occurrences_reelles = 15;
    occurrences_aleatoires = 15;
    fold_change = calculer_fold_change(occurrences_reelles, occurrences_aleatoires);
    printf("Fold change pour %d occurrences réelles et %d occurrences aléatoires : %.2f\n", occurrences_reelles, occurrences_aleatoires, fold_change);
    assert(fold_change > 0.99 && fold_change < 1.01);  // Fold change ≈ 1.0
    
    occurrences_reelles = 5;
    occurrences_aleatoires = 10;
    fold_change = calculer_fold_change(occurrences_reelles, occurrences_aleatoires);
    printf("Fold change pour %d occurrences réelles et %d occurrences aléatoires : %.2f\n", occurrences_reelles, occurrences_aleatoires, fold_change);
    assert(fold_change > 0.49 && fold_change < 0.51);  // Fold change ≈ 0.5
    
    occurrences_reelles = 10;
    occurrences_aleatoires = 0;
    fold_change = calculer_fold_change(occurrences_reelles, occurrences_aleatoires);
    printf("Fold change pour %d occurrences réelles et %d occurrences aléatoires : %.2f\n", occurrences_reelles, occurrences_aleatoires, fold_change);
    assert(fold_change > 999999.0);  // Fold change très grand à cause de l'epsilon
    
    printf("Test de calculer_fold_change passé avec succès.\n");
}

// Test programme principale
void test_etendre_et_traiter_k_uplets() {
    // Séquence génomique de test (courte pour la visualisation)
    char *sequence_reelle = "ATGCCTGCACTGCGATCGT";  // Séquence réelle simulée
    char *sequence_alea = "CGTACGATCGTACGTACGTA";   // Séquence aléatoire simulée (pour comparer)
    int longueur_sequence_reelle = strlen(sequence_reelle);
    
    // Simuler la position d'un gène dans la séquence réelle
    int position_gene = 15;  // Position du gène simulée dans la séquence réelle

    // Vérifier que la région en amont ne dépasse pas les limites de la séquence
    int taille_region = 5;  // Taille de la région en amont à analyser
    if (position_gene - taille_region < 0) {
        printf("Erreur : la région en amont dépasse les limites de la séquence\n");
        return;
    }

    // Extraire la sous-séquence à analyser en amont du gène
    char* sous_sequence = (char*)malloc((taille_region + 1) * sizeof(char));
    if (sous_sequence == NULL) {
        printf("Erreur d'allocation mémoire pour la sous-séquence\n");
        return;
    }
    strncpy(sous_sequence, &sequence_reelle[position_gene - taille_region], taille_region);
    sous_sequence[taille_region] = '\0';

    // Paramètres pour l'analyse des k-uplets
    int longueur_k_uplet = 2;     // Longueur des k-uplets à rechercher
    int longueur_min_motif = 4; // Longueur minimale d'un motif pour le test statistique

    // Appeler la fonction pour traiter les k-uplets
    printf("=== Test 1 : Appel de traiter_k_uplets avec des séquences courtes ===\n");
    traiter_k_uplets(sous_sequence, taille_region, sequence_reelle,
                     longueur_sequence_reelle, longueur_k_uplet, longueur_min_motif,
                     sequence_alea, position_gene);

    // Test 2 : Extension d'un motif particulier après analyse
    printf("\n=== Test 2 : Test de l'extension d'un motif ===\n");
    Motif motif_test;
    motif_test.sequence = strdup("TGCA");
    motif_test.start_pos = 4;  // Position du motif "TGCA" dans la séquence
    motif_test.length = strlen(motif_test.sequence);
    int occurrences_reelles = 1;  // Simuler une occurrence

    int extension_result = etendre_k_uplet(sequence_reelle, longueur_sequence_reelle, &motif_test, &occurrences_reelles);
    printf("Résultat de l'extension : %d\n", extension_result);
    printf("Motif après extension : %s\n", motif_test.sequence);
    printf("Position du motif après extension : %d\n", motif_test.start_pos);
    printf("Longueur du motif après extension : %d\n", motif_test.length);
    printf("Occurrences réelles après extension : %d\n", occurrences_reelles);

    // Assertions pour vérifier le comportement attendu
    assert(strcmp(motif_test.sequence, "CTGCA") == 0);  // Test attendu après l'extension
    assert(motif_test.length == 5);

    free(motif_test.sequence);

    // Libérer la mémoire
    free(sous_sequence);

    printf("\n=== Fin des tests combinés de etendre_k_uplet et traiter_k_uplets ===\n");
}