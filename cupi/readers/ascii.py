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

def read_cyc(fcyc):
  with open(fcyc, 'r') as f:
    lines = f.readlines()
  # read particle labels
  labels = []
  for iline, line in enumerate(lines):
    tokens = line.split()
    if tokens[0] == 'type':
      labels.append(tokens[1])
    else:
      break
  # read cycle lengths
  data = []
  for line in lines[iline:]:
    tokens = line.split()
    entry = {'ncycle': int(tokens[0])}
    for itype, label in enumerate(labels):
      entry.update({label: float(tokens[itype+1])})
    data.append(entry)
  return data
