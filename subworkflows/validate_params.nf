process validate_params {

    input:
    val pathsToValidate
    
    script:
    """

    IFS=',' read -ra paths <<< "${pathsToValidate}"

    for path in "\${paths[@]}"; do
        if [ -e "\${path}" ]; then
            continue
         else
             echo "Error: '\${path}' does not exist" 
             exit 1
        fi
    done
    
    """
}