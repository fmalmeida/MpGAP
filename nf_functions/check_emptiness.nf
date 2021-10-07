def is_empty(in_ch) {
    file('my_file.txt', checkIfExists: false).write('ALL OK\n')
    in_ch.ifEmpty { file('my_file.txt', checkIfExists: false).write('EMPTY\n') }
    result    = file('my_file.txt', checkIfExists: false).readLines().findAll { it.contains('EMPTY') }.toString()
    if (result == '[EMPTY]') {
        empty = 'true'
      } else {
        empty = 'false'
      }
    file('my_file.txt', checkIfExists: false).delete()
    return empty
}