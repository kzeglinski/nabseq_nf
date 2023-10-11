process check_report_group {
    input:
    path sample_sheet

    script:
    """
    #!/usr/bin/env python

    import csv

    with open("${sample_sheet}", 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        if 'report_group' not in header:
            # Replace the header with the new header
            new_header = header + ['report_group']
            rows = [new_header]
            for row in reader:
                row.append('1')
                rows.append(row)

            # Write the modified data back to the sample sheet
            with open("${sample_sheet}", 'w') as f:
                writer = csv.writer(f)
                writer.writerows(rows)

    """
}

process check_whitespace {
    input:
    path sample_sheet

    script:
    """
    # Check whitespace existence

    if grep -q ' ' "${sample_sheet}"; then
        echo "Error: Whitespace detected in the sample sheet. Remove before continuing."
        exit 1
    fi

    """
}

process check_header {
    input:
    path sample_sheet

    script:
    """
     # Check Header 

    expected_header="barcode,sample_name,species,report_group"
    header=\$(head -n 1 "${sample_sheet}" | sed 's/.\$//')

    if [ "\${header}" != "\${expected_header}" ]; then
        echo "Error: 
              Expected Header: barcode,sample_name,species,report_group
              Current Header : \${header}"
        exit 1
    fi
    """
}

process check_reference_sequences {
    input:
    path reference_dir

    script:
    """
    if [ -z "\$(ls -A "${reference_dir}")" ]; then
        echo "Error: The 'reference_sequences' directory is empty."
        exit 1
    fi
    """
}


process check_references {
    input:
    path sample_sheet
    path reference_dir

    script:
    """
    #!/bin/bash

    species_list=()
    first_line=true
    missing_species=""

    while IFS=',' read -r _ _ species _; do
        # Skip the first line (header)
        if [[ \$first_line == true ]]; then
            first_line=false
            continue
        fi

        # Check if the species value is not empty and not already in the list
        if [[ -n "\$species" && ! "\${species_list[@]}" =~ "\$species" ]]; then
            species_list+=("\$species")
        fi
    done < "${sample_sheet}"

    reference_list=()

    # Loop through files in the directory
    for file in "${reference_dir}"/*_all_imgt_refs.fasta; do
        # Extract the animal name from the filename (first word)
        animal=\$(basename "\$file" | cut -d'_' -f1)
        
        # Check if the animal is not already in the list
        if [[ ! " \${reference_list[@]} " =~ " \${animal} " ]]; then
            reference_list+=("\$animal")
        fi
    done

    for species in "\${species_list[@]}"; do
        if [[ ! " \${reference_list[@]} " =~ " \${species} " ]]; then
            if [ -z "\$missing_species" ]; then
                missing_species="\$species"
            else
                missing_species="\${missing_species}, \$species"
            fi
        fi
    done

    if [ -n "\$missing_species" ]; then
        echo "Error: Species listed in the sample sheet do not have reference sequences in the reference directory. 
        Missing species: \$missing_species"
        exit 1
    fi

    """
}


workflow validate_sample_sheet {

    take:
        sample_sheet
        reference_dir

    main:
        check_report_group(sample_sheet)
        check_whitespace(sample_sheet)
        check_header(sample_sheet)
        check_reference_sequences(reference_dir)
        check_references(sample_sheet, reference_dir)
}