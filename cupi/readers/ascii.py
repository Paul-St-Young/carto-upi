def read_acc(facc):
  with open(facc, 'r') as f:
    lines = f.readlines()
  assert ('gate' in lines[0])
  assert ('nsets' in lines[-1])

  data = []
  for line in lines[1:-1]:
    tokens = line.split()
    # first collect columns that are always there
    labels = ['ptype', 'mvtype', 'nmover', 'racc', 'ntry', 'nacc']
    dtypes = [str, str, int, float, int, int]
    ncol = len(labels)
    entry = dict()
    for val, label, dtype in zip(tokens, labels, dtypes):
      entry[label] = dtype(val)
    # next look at what is left
    extra_tokens = tokens[ncol:]
    nleft = len(extra_tokens)
    if nleft > 0:  # should be level acceptance
      lacc = map(float, extra_tokens)
      for ilvl, racc in enumerate(lacc):
        name = 'l%d_racc' % ilvl
        entry.update({name: racc})
    data.append(entry)
  return data
