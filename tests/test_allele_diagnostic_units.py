from chronoclade.temporal_diagnostics import date_regression, allele_unit_diagnostic


def test_allele_unit_scale_keeps_fraction_evidence_and_residuals():
    diagnostic = date_regression([dict(sample_id=str(i), collection_date=str(2000+i), d=i*.01)
                                  for i in range(3)], distance_field='d',
                                 distance_units='fraction', label='test')
    scaled = allele_unit_diagnostic(diagnostic, 100)
    assert scaled['slope'] == diagnostic['slope'] * 100
    assert scaled['fraction_slope'] == diagnostic['slope']
    assert scaled['points'][1]['d'] == 1
    assert diagnostic['points'][1]['d'] == .01
    assert 'scheme-equivalent' in scaled['slope_units']
    fixed = allele_unit_diagnostic(diagnostic, 100, fixed_callable=True)
    assert 'fixed callable' in fixed['distance_units']
