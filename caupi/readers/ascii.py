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
  def get(tag, ival, dtype, pick=None):
    for line in lines:
      if line.strip().startswith(tag):
        tokens = line.split()
        if pick is not None:
          ipos, myval = pick
          if tokens[ipos] != myval:
            continue
        val_text = tokens[ival]
    try:
      val = dtype(val_text)
    except UnboundLocalError:
      msg = '%s not in %s' % (tag, text)
      raise RuntimeError(msg)
    except ValueError:
      print(tag, text)
    return val
  # read temperature
  beta = get('BETA', 1, float)
  nslice = get('NSLICES', 1, int)
  # read electrons
  lam = get('TYPE', 2, float, pick=(1, 'e'))
  itype = get('TYPE', 3, int, pick=(1, 'e'))
  nelec = get('TYPE', 4, int, pick=(1, 'e'))
  entry = {'beta': beta, 'nslice': nslice, 'nelec': nelec, 'lam': lam}
  if itype > 1:
    nup = get('TYPE', 5, int, pick=(1, 'e'))
    ndn = get('TYPE', 6, int, pick=(1, 'e'))
    assert nup+ndn == nelec
    entry.update({'nup': nup, 'ndn': ndn})
  # read box
  lx = get('BOXSIZE', 1, float)
  ly = get('BOXSIZE', 2, float)
  lz = get('BOXSIZE', 3, float)
  entry['box'] = (lx, ly, lz)
  # read pairpot
  kcut = get('CUTK', 1, float)
  ifewald = get('CUTK', 2, int)
  entry['kcut'] = kcut
  entry['ifewald'] = ifewald
  # read fixnode
  try:
    nref = get('NREF', 1, int)
  except RuntimeError as err:
    if str(err).startswith('NREF not in '):
      nref = 0
  try:
    ifn = get('FIXNODE', 1, int)
  except RuntimeError as err:
    if str(err).startswith('FIXNODE'):
      ifn = 0
    else:
      raise err
  try:
    get('NOLEAK', 0, str)
    ifleak = 1
  except RuntimeError as err:
    if str(err).startswith('NOLEAK'):
      ifleak = 0
    else:
      raise err
  entry['nref'] = nref
  entry['ifn'] = ifn
  entry['ifleak'] = ifleak
  # anything I missed
  entry['text'] = text
  return entry

def read_sy(fsy):
  with open(fsy, 'r') as f:
    text = f.read()
  return parse_sy(text)
