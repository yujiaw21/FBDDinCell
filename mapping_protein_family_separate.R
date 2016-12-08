target.list <- read.csv("target.list.csv", header = TRUE)

terms <-
  c(
    "adatptors",
    "channels",
    "chaperones",
    "dna_rna_bind",
    "enzymes",
    "fatty_acid_lip_meta",
    "heme_related",
    "lipid_bind",
    "nucbinding",
    "receptors",
    "repressors",
    "structural_proteins",
    "tf_and_regulators",
    "transporter",
    "no_category"
  )

KW_no_category <-
  read.delim(
    "./protein_families/families/no_category.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )
KW_adatptors <-
  read.delim(
    "./protein_families/families/KW_adatptors.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )
KW_channels <-
  read.delim(
    "./protein_families/families/KW_channels.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )
KW_chaperones <-
  read.delim(
    "./protein_families/families/KW_chaperones.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )
KW_dna_rna_bind <-
  read.delim(
    "./protein_families/families/KW_dna_rna_bind.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )
KW_enzymes <-
  read.delim(
    "./protein_families/families/KW_enzymes.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )
KW_fatty_acid_lip_meta <-
  read.delim(
    "./protein_families/families/KW_fatty_acid_lip_meta.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )
KW_heme_related <-
  read.delim(
    "./protein_families/families/KW_heme_related.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )
KW_lipid_bind <-
  read.delim(
    "./protein_families/families/KW_lipid_bind.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )
KW_nucbinding <-
  read.delim(
    "./protein_families/families/KW_nucbinding.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )
KW_receptors <-
  read.delim(
    "./protein_families/families/KW_receptors.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )
KW_repressors <-
  read.delim(
    "./protein_families/families/KW_repressors.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )
KW_structural_proteins <-
  read.delim2(
    "./protein_families/families/KW_structural_proteins.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )
KW_tf_and_regulators <-
  read.delim(
    "./protein_families/families/KW_tf_and_regulators.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )
KW_transporter <-
  read.delim(
    "./protein_families/families/KW_transporter.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )

for (i in 1:length(terms)) {
  target.list[, terms[i]] <- ""
  target.list[target.list$Identifier %in% eval(parse(text = paste("KW_", terms[i], sep =
                                                                    "")))$V1, terms[i]] <- "*"
}



#target.list[target.list$Identifier %in% KW_adatptors$V1,]$adatptors <- "*"

write.csv(target.list, file = "target.list.w.protein_families.csv", row.names = FALSE)
