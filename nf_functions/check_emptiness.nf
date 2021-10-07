def is_empty(in_ch, variable) {
    in_file = file("${variable}.txt", checkIfExists: false)
    in_ch.ifEmpty { in_file.write('EMPTY\n') }
    result = in_file.readLines().findAll { it.contains('EMPTY') }.toString()
    if (result == '[EMPTY]') {
        empty = 'true'
      } else {
        empty = 'false'
      }
    file("${variable}.txt", checkIfExists: false).delete()
    return empty
}