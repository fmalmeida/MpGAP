def is_empty(in_ch, variable) {
    // the if empty conditional does not create the file in the global pattern
    // thus, when checking if the file exists, it will return true only when channel is loaded
    in_file = file("${variable}.txt", checkIfExists: false)
    in_ch.ifEmpty { in_file.write('EMPTY') }
    if (in_file.exists()) {
      empty = 'OK'
    } else {
      empty = 'EMPTY'
    }
    return empty
}