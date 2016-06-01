from hindsight.bibtex.bib import split_on_comma_or_space

def test_split_on_comma_or_space():
    assert set(split_on_comma_or_space('pflB adhE frdA')) == set(['pflB', 'adhE', 'frdA'])
    assert set(split_on_comma_or_space('pflB, adhE,  frdA')) == set(['pflB', 'adhE', 'frdA'])
    assert set(split_on_comma_or_space('pflB, adhE,  frdA, ')) == set(['pflB', 'adhE', 'frdA'])
