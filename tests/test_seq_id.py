#!/usr/bin/env python3
'''
Calculate the sequence identity between the sequences in an alignment
'''

import seq_id


def test_calc_seqids():
    '''
    Test calc_dist()
    '''
    sequences = ['AAAAACCCCC', 'AAAAAAAAAA',
                 'AACCCCCCCC', 'CCCCCAAAAA']
    expected = [0.5, 0.7, 0, 0.2, 0.5, 0.3]
    assert seq_id.calc_seqids(sequences) == expected


def test_compare_res():
    '''
    Test compare_res()
    '''
    # equal sequences
    seq1 = 'ACDEFGHIKLMNPQRSTVWY'
    seq2 = 'ACDEFGHIKLMNPQRSTVWY'
    expected = [1]*20
    assert seq_id.compare_res(seq1, seq2) == expected
    # different sequences
    seq1 = 'AAAAAAAAAAAAAAAAAAAY'
    seq2 = 'ACDEFGHIKLMNPQRSTVWY'
    expected = [1] + [0]*18 + [1]
    assert seq_id.compare_res(seq1, seq2) == expected
    # skip gaps
    seq1 = 'ACDEF-----MNPQRSTVWY'
    seq2 = 'ACDEFGHIKLMNPQRSTVWY'
    expected = [1]*15
    assert seq_id.compare_res(seq1, seq2) == expected
    # gaps = missmatch
    expected = [1]*5 + [0]*5 + [1]*10
    assert seq_id.compare_res(seq1, seq2, False) == expected


def test_main(capsys):
    '''
    Test main()
    '''
    import tempfile

    indata = ('>seq1\n' +
              'MLVSRRLTGARARAPLLASLLEAWCRQGRTTSSYSAFSEPSHVRALVYGNHGDPAKVIQLKNLELTAVEG' +
              'SDVHVKMLAAPINPSDINMIQGNYGLLPKLPAVGGNEGVGQVIAVGSSVSGLKPGDWVIPANAGLGTWRT' +
              'EAVFSEEALIGVPKDIPLQSAATLGVNPCTAYRMLVDFEQLQPGDSVIQNASNSGVGQAVIQIASALGLK' +
              'GGTMVTYGGMSKKPITVSTTSFIFKDLALRGFWLQSWLSMGKVKECREMIDYLLGLARDGKLKYETELVP\n' +
              '>seq2\n'
              'KENDVCVKMIAAPINPSDINRIEGVYPVRPPVPAVGGYEGVGEVYAVGSNVNGFSPGDWVIPSPPSSGTW' +
              'QTYVVKEESVWHKIDKECPMEYAATITVNPLTALRMLEDFVNLNSGDSVVQNGATSIVGQCVIQLARLRG' +
              'ISTINLIRDRAGSDEAREQLKALGADEVFSESQLNVKNVKSLLGNLPEPALGFNCVGGNAASLVLKYLRE' +
              'GGTMVTYGGMSKKPITVSTTSFIFKDLALRGFWLQSWLSMGKVKECREMIDYLLGLARDGKLKYETELVP\n' +
              '>seq3\n'
              'KENDVCVKMIAAPINPSDINRIEGVYPVRPPVPAVGGYEGVGEVYAVGSNVNGFSPGDWVIPSPPSSGTW' +
              'QTYVVKEESVWHKIDKECPMEYAATITVNPLTALRMLEDFVNLNSGDSVVQNGATSIVGQCVIQLARLRG' +
              'EAVFSEEALIGVPKDIPLQSAATLGVNPCTAYRMLVDFEQLQPGDSVIQNASNSGVGQAVIQIASALGLK' +
              'GGTMVTYGGMSKKPITVSTTSFIFKDLALRGFWL-------------EMIDYLLGLARDGKLKYETELVP\n')

    fasta_filename = tempfile.mkstemp()[1]
    with open(fasta_filename, 'w') as tmpf:
        tmpf.write(indata)

    # test normal run
    assert seq_id.main(fasta_filename) is None
    out, err = capsys.readouterr()
    expected = 'Median: 0.5243\nAverage: 0.5342\nMax: 0.7603\nMin: 0.3179\n'
    assert out == (expected)
    # test incorrect config
    assert seq_id.main(fasta_filename, ['qq']) == 'Unknown option: qq'
    assert seq_id.main(fasta_filename, ['gm']) is None
    out, err = capsys.readouterr()
    expected = 'Median: 0.5\nAverage: 0.5143\nMax: 0.725\nMin: 0.3179\n'
    assert out == (expected)
    assert seq_id.main(fasta_filename, ['sn']) is None
    out, err = capsys.readouterr()
    assert out == ''


def test_print_stats(capsys):
    '''
    Test print_stats()
    '''
    data = [0.3, 0.65, 0.2, 0.75, 0.6]
    seq_id.print_stats(data)
    out, err = capsys.readouterr()
    assert out == ('Median: 0.6\nAverage: 0.5\n' +
                   'Max: 0.75\nMin: 0.2\n')
    # test 4 decimal digits
    data = [0.1234, 0.1234, 0.12345, 0.12345, 0.54321]
    seq_id.print_stats(data)
    out, err = capsys.readouterr()
    assert out == ('Median: 0.1235\nAverage: 0.2074\n' +
                   'Max: 0.5432\nMin: 0.1234\n')


def test_set_config():
    '''
    Test set_config()
    '''
    # one option
    assert seq_id.set_config(['gm']) == [False, True]
    # incorrect option
    assert seq_id.set_config(['as']) == 'Unknown option: as'
    # multiple options, "incorrect" order
    assert seq_id.set_config(['sn', 'gs']) == [True, False]
    # should use the last given parameter if multiple copies are given
    assert seq_id.set_config(['gs', 'gs', 'gm']) == [False, True]
