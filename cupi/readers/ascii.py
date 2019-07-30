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

def parse_sy(text):
  lines = text.split('\n')
  def get(tag, ival, dtype):
    for line in lines:
      if line.strip().startswith(tag):
        val_text = line.split()[ival]
    try:
      val = dtype(val_text)
    except UnboundLocalError:
      msg = '%s not in %s' % (tag, text)
      raise RuntimeError(msg)
    except ValueError:
      print(tag, text)
    return val
  beta = get('BETA', 1, float)
  nslice = get('NSLICES', 1, int)
  lam = get('TYPE', 2, float)
  itype = get('TYPE', 3, int)
  nelec = get('TYPE', 4, int)
  try:
    nref = get('NREF', 1, int)
  except RuntimeError as err:
    if str(err).startswith('NREF not in '):
      nref = 0
  entry = {'beta': beta, 'nslice': nslice, 'nref': nref, 'nelec': nelec, 'lam': lam}
  if itype > 1:
    nup = get('TYPE', 5, int)
    ndn = get('TYPE', 6, int)
    assert nup+ndn == nelec
    entry.update({'nup': nup, 'ndn': ndn})
  return entry

def read_sy(fsy):
  with open(fsy, 'r') as f:
    text = f.read()
  return parse_sy(text)
