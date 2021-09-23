process medaka_models {
    output:
    stdout models
    label 'main'

    """
    source activate MEDAKA;
    medaka tools list\_models | grep "Available"
    """

}