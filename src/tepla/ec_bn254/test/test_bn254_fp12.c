#include <assert.h>

#include "rdtsc.h"

#include "../ec_bn254_lcl.h"

#define N 100
#define M 10

//============================================
//   Feature of Field
//============================================
void test_feature(Field f)
{
    fprintf(stdout, "---\n");
    fprintf(stdout, "FieldType: %s\n", field_get_name(f));
    gmp_fprintf(stdout, "   characteristic: %Zx\n", *field_get_char(f));
    gmp_fprintf(stdout, "   order of field: %Zx^%d\n", *field_get_char(f), field_get_degree(f));
    fprintf(stdout, "---\n");
}

//============================================
//   test for arithmetic operations
//============================================
void test_arithmetic_operation_beuchat(Field f)
{
    int i;
    unsigned long long int t1, t2;
    Element a, b, c, d;
    Element d1, d2;
    Element h1, h2, h3;
    Element g1, g2, g3;

    mpz_t exp;

    //--------------------
    //  init
    //--------------------
    element_init(a, f);
    element_init(b, f);
    element_init(c, f);
    element_init(d, f);

    element_init(d1, f->base);
    element_init(d2, f->base);

    element_init(h1, f->base->base);
    element_init(h2, f->base->base);
    element_init(h3, f->base->base);

    element_init(g1, f->base->base);
    element_init(g2, f->base->base);
    element_init(g3, f->base->base);

    //--------------------
    //  add
    //--------------------
    element_set_str(a, "1F8989298CBA4F1D8530E545C86D5AC9C1E04EA514305C4A5742E156B7C39CD1 1D1FEC91A6D00B497F846D80202AFC27184141E26051C63C3EC8974721E1EA78 7ED4FD8F56E3F1C3D9BB1832F4A66A5CCEE2976851FF781E7656FCA26D42BCE E8D93F6E0661803549F4CE03854FF513BD900005D97FEF90A38B1D28260FFDC 156676170BA7ECB549209DD447C3A06DA8AC548057399CA5E1730FF7380C2D48 57516260F5BCB09714F2032E615C24C029D87F966524F46A833231649FD0CBE 1123380E65F1BB40258BD8CC829C2953290D4D54191A0B7DA0075DE10CF48A80 8DCDAECFD365670E30F416A957060E51D94766D8AD951450917A08A88402273 2017BBE296901940107189E623C0CEA14B3DBDF068D6CF7D643A49E314E48AF0 B21605AD67C6026A4AD7D5A4DCE66FF2D3A7A15FA3E45234A17A76D5C4A5624 1090186BC36CA99F9C04D613E10693CA93C058CF8A20AED0CFC138C291B17C68 19399C1D77C66053327286ADC5851EC4C18E5008F4D7C2B02818EAC2B29A99B");
    element_set_str(b, "101956683C21E941221B1440219BEAD08EC535816EB2A73FA45E9F2268909B41 EAEEC005EB27B3FD4B141D6F40C75D30E216CCE00CF25C8E6AA0520C08ED74F 21C4065DBDD1EC0FBB6E4A0FC38F815D0A5593454D9885C065BC93A2927D1C5D 1A938F7BF4CF6D6B770F807C64C1663252DE21FB3495094A9EBB8EA942871042 DFACA46B7E39E919FBCC652BD0DAEBEDCE2ABFF3BA6ACC91CE877099AE71C35 1679AEFC82A0903AB73AD7154E262106F05AF2708D708ECCE7DCBF39EB35D4FE 1E907D4E9CB0EAEF8CCC438830FA28EEDFCD5190A828F21812E58D2AD526DF0A 135DEF8CB2DDFE5A413F5FDAFC52E2FDEE5D0AC215B44226BFEA50905E5CC03B 94040DF603958725183122D4CE4BC052FD5CBFC9019ECF7CFD5DC1429A8C39B 2112DCA6D5CD9E4403D7C4E2098761479EC3D209FE6AC985C8FEAF8D0FC3B658 17E7D50EE673800EBCF9E91AA7F55FECBEBC43AAE4451535E35082A95252DF2B 1E1CAA720242AD9C939194FF68DD65A674F17D8DE06BD9F0C7D5F3844B602C64");
    element_set_str(d, "C31E48D2B9B28A058D5DEED63244358D3A344E66AE303887D21207920543811 85DDD8D684176CB05BF94BE8D526FB8A9606F704920EC03A6F23C67E270C1C6 6405B3215FF1B6DAA93E0FA6BF4E5C15A417D7BBAB87D40CEA1A36CB951482A 5B0286E37F475B07D38B2C41631634211B4E2BB7A2D08422A73E07BC4E8101D 2361405DC38B8B46E8DD642704D14F2C858F007F92E0496EFE5B8700D2F3497D 1BEEC52291FC5B442889F748343BE352F2F87A69F3C2DE13900FE2503532E1BC C42BA586561967163E201BC2CB150008BD85FA4A942FD94346C8B0BE21B6989 1C3ACA79B01454CB244EA14591C343E30BF1812FA08D936BC901F11AE69CE2AE 5E701BD598861F4137E817AE9C08864FE114AACE0F0BC73B58FC5F73E8D4E8A 8C341FD0F08EEAC5A0F27A3D070C6054EFC0CDFE0A90EA79495F6FA6C0E0C7B 506F2760C9F19F00A88A4960216F175D57A5D3A5665C40534915B6BE4045B92 1FB04433D9BF13A1C6B8BD6A4535B792C10A628E6FB9561BCA5782307689D5FF");

    element_add(c, a, b);

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_add(c, a, b);
    }
    t2 = rdtsc();

    printf("element add: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //--------------------
    //  sub
    //--------------------
    element_set(d, c);
    element_sub(c, c, d);

    assert(element_is_zero(c));

    //--------------------
    //  mul
    //--------------------
    element_mul(c, a, b);

    element_set_str(h1, "D1A561303CE3F6C279C549E37F06204EE6FC5DB520CB879EF3F0D69732E1FFD 1CBD9E049955F69C23733968E50CEA6550CA95AE91943DF417651F01E2610F11");
    element_set_str(h2, "1035E673CFA033AC18614540CC49E6A6F62BFC1DFE5FED8A346C5A562C359AC5 2201664E07EAF1F77C63BC20469FB1B8610CA809DB0EFC7AE3E0E5CDFBC4631A");
    element_set_str(h3, "2073EE2024349E3B9EF988F59CC71AC3A0D1C9CCC4E47E42129DEFB6BE42C3D 113490FF65772E0A7C04D284AA35B44D7D46AA0AB54F81C570EB13C241850219");

    element_set_str(g1, "2FE2518500A64C7E98EAC34A757FABC7F216168212051BB480D99400C8C41B2 1F2462AF9BC4E3C3C2046345265B175DAE717CF51CE3914AA99CA3BCDE33DAF1");
    element_set_str(g2, "1D20A00EE8D645643F2AA82F89C6E57FC37D24ACC5A03BB8D470F4536D3ED98F BEEA0DAD730708B89B6E80521CC36000D74697D720ECC637DF999E14D447CD1");
    element_set_str(g3, "91EB3EC79C0EDFAAEDE8BC17B5737F6C84A26A3E29A4E68B24D74862D67CC9C 1AAD97B2CEB4A493416FD27B164EAF20D4A32E28574C50857E2EFC7D11D18CE0");

    element_set(((Element *)d1->data)[0], h1);
    element_set(((Element *)d1->data)[1], h2);
    element_set(((Element *)d1->data)[2], h3);

    element_set(((Element *)d2->data)[0], g1);
    element_set(((Element *)d2->data)[1], g2);
    element_set(((Element *)d2->data)[2], g3);

    element_set(((Element *)d->data)[0], d1);
    element_set(((Element *)d->data)[1], d2);

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_mul(c, a, b);
    }
    t2 = rdtsc();

    printf("element mul: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //--------------------
    //  sqr
    //--------------------
    element_sqr(c, a);
    element_mul(d, a, a);

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_sqr(c, a);
    }
    t2 = rdtsc();

    printf("element sqr: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //--------------------
    //  random
    //--------------------
    element_random(b);

    //--------------------
    //  inv
    //--------------------
    element_mul(c, a, b);
    element_inv(b, b);
    element_mul(c, c, b);
    element_inv(d, a);
    element_mul(d, a, d);

    assert(element_cmp(c, a) == 0);
    assert(element_is_one(d));

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_inv(c, a);
    }
    t2 = rdtsc();

    printf("element inv: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //--------------------
    //  pow
    //--------------------
    mpz_init(exp);

    mpz_set_str(exp, "103DC9103A8BDD2D08BA8E58EC4AD8AA108177A04CFBBAE6FB7CE07A88E298C3", 16);

    element_pow(c, a, exp);
    element_set_str(d, "7069546767E6FD9CD909387BCFBB30BA724199B62FBD76937B793E782516E4C 2923EA6D063196EE83E9D381BDF49EEBE9752738059298078DD054AE91C598F 110AA75DB2914DAC5A0FBBCE56719CD5F95D40FF65744E4307C3EEDDBDF02453 5F4ABCA5A37579FF1EE2E07AD2201C61E4E10D5B95FA4BE3A03D98ED7EB3F10 124955B6581D6A93DACBB0982598475F17E3E718AE98FB0EA859547022B57AA5 1A78997B1C862DC05C742EEEDA0FA036351306CE2CF168DC9EBF8DCDC5814BA2 B930A73F8505797959C1DDCA947AABD4FC356870DBD0CDE87301C50E40AFD15 140FB498FBA6110D1A6D13CF7D09400FDC50AB73CC8A30AAE8F7893A550F9404 82ABDE131431BB692F4FD80D7BB5F90C8AF4B339396A738FB124E4295C37002 1D37E5244B85179C4DB05C7160DB55DD98661FF69BB6F090969EFDDEABDB9A28 10819CF9418EA94B13FD687CB1E1E73AAF478BC2BFA5474E50CD17AE63EB6618 EA8DDA68441555D8367BBBBA5BB3719E2C0849BD4CE0D95E262F0D22719587A");

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_pow(b, a, exp);
    }
    t2 = rdtsc();

    printf("element pow with torsion: %.2lf [clock]\n", (double)(t2 - t1) / N);

    mpz_set(exp, f->order);

    for (i = 0; i < 10; i++)
    {
        element_random(a);

        element_pow(b, a, exp);

        assert(element_cmp(b, a) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        element_pow(b, a, exp);
    }
    t2 = rdtsc();

    printf("element pow with order: %.2lf [clock]\n", (double)(t2 - t1) / M);

    mpz_clear(exp);

    //--------------------
    //  clear
    //--------------------
    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(d);

    element_clear(d1);
    element_clear(d2);

    element_clear(h1);
    element_clear(h2);
    element_clear(h3);
    element_clear(g1);
    element_clear(g2);
    element_clear(g3);
}

void test_arithmetic_operation_aranha(Field f)
{
    int i;
    unsigned long long int t1, t2;
    Element a, b, c, d;
    Element d1, d2;
    Element h1, h2, h3;
    Element g1, g2, g3;

    mpz_t exp;

    //--------------------
    //  init
    //--------------------
    element_init(a, f);
    element_init(b, f);
    element_init(c, f);
    element_init(d, f);

    element_init(d1, f->base);
    element_init(d2, f->base);

    element_init(h1, f->base->base);
    element_init(h2, f->base->base);
    element_init(h3, f->base->base);

    element_init(g1, f->base->base);
    element_init(g2, f->base->base);
    element_init(g3, f->base->base);

    //--------------------
    //  add
    //--------------------
    element_set_str(a, "1F8989298CBA4F1D8530E545C86D5AC9C1E04EA514305C4A5742E156B7C39CD1 1D1FEC91A6D00B497F846D80202AFC27184141E26051C63C3EC8974721E1EA78 7ED4FD8F56E3F1C3D9BB1832F4A66A5CCEE2976851FF781E7656FCA26D42BCE E8D93F6E0661803549F4CE03854FF513BD900005D97FEF90A38B1D28260FFDC 156676170BA7ECB549209DD447C3A06DA8AC548057399CA5E1730FF7380C2D48 57516260F5BCB09714F2032E615C24C029D87F966524F46A833231649FD0CBE 1123380E65F1BB40258BD8CC829C2953290D4D54191A0B7DA0075DE10CF48A80 8DCDAECFD365670E30F416A957060E51D94766D8AD951450917A08A88402273 2017BBE296901940107189E623C0CEA14B3DBDF068D6CF7D643A49E314E48AF0 B21605AD67C6026A4AD7D5A4DCE66FF2D3A7A15FA3E45234A17A76D5C4A5624 1090186BC36CA99F9C04D613E10693CA93C058CF8A20AED0CFC138C291B17C68 19399C1D77C66053327286ADC5851EC4C18E5008F4D7C2B02818EAC2B29A99B");
    element_set_str(b, "101956683C21E941221B1440219BEAD08EC535816EB2A73FA45E9F2268909B41 EAEEC005EB27B3FD4B141D6F40C75D30E216CCE00CF25C8E6AA0520C08ED74F 21C4065DBDD1EC0FBB6E4A0FC38F815D0A5593454D9885C065BC93A2927D1C5D 1A938F7BF4CF6D6B770F807C64C1663252DE21FB3495094A9EBB8EA942871042 DFACA46B7E39E919FBCC652BD0DAEBEDCE2ABFF3BA6ACC91CE877099AE71C35 1679AEFC82A0903AB73AD7154E262106F05AF2708D708ECCE7DCBF39EB35D4FE 1E907D4E9CB0EAEF8CCC438830FA28EEDFCD5190A828F21812E58D2AD526DF0A 135DEF8CB2DDFE5A413F5FDAFC52E2FDEE5D0AC215B44226BFEA50905E5CC03B 94040DF603958725183122D4CE4BC052FD5CBFC9019ECF7CFD5DC1429A8C39B 2112DCA6D5CD9E4403D7C4E2098761479EC3D209FE6AC985C8FEAF8D0FC3B658 17E7D50EE673800EBCF9E91AA7F55FECBEBC43AAE4451535E35082A95252DF2B 1E1CAA720242AD9C939194FF68DD65A674F17D8DE06BD9F0C7D5F3844B602C64");
    element_set_str(d, "A7F7B0F88DC385CED17AC05EA094591EF84842682E3037654A18079205437FF 6AB740FC58286879A0161D7143771F1C541AEB06120EBF17E729C67E270C1B4 48DF1B473402B2A3ED5AE12F2D9E7FA7622BCBBD2B87D2EA622036CB9514818 3FDBEF09535856D117A7FDC9D16657B2D9621FB922D083001F4407BC4E8100B 2361405DC38B8B46E8DD642704D14F2C858F007F92E0496EFE5B8700D2F3497D 1BEEC52291FC5B442889F748343BE352F2F87A69F3C2DE13900FE2503532E1BC A9050DAC2A2A62DF823CED4B3965239A7B99EE4C142FD820BECEB0BE21B6977 1C3ACA79B01454CB244EA14591C343E30BF1812FA08D936BC901F11AE69CE2AE 434983FB6C971B0A7C04E9370A58A9E19F289ECF8F0BC618D1025F73E8D4E78 710D87F6C49FE68EE50F4BC5755C83E6ADD4C1FF8A90E956C1656FA6C0E0C69 35488F869E029AC9ECA71AE88FBF3AEF15B9C7A6E65C3F30C11BB6BE4045B80 1FB04433D9BF13A1C6B8BD6A4535B792C10A628E6FB9561BCA5782307689D5FF");

    element_add(c, a, b);

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_add(c, a, b);
    }
    t2 = rdtsc();

    printf("element add: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //--------------------
    //  sub
    //--------------------
    element_set(d, c);
    element_sub(c, c, d);

    assert(element_is_zero(c));

    //--------------------
    //  mul
    //--------------------
    element_mul(c, a, b);

    element_set_str(h1, "A1468F71AEF6E895F1BA2C43F6167137EDC410B7E939D2EC4257344582ADAC9 19E5A451176844970B07F64CFE4E638475295458B3F54E746FD5B3EE41C13C08");
    element_set_str(h2, "12B4050B9DEEDCFCB298B3E4A12E4C2FD10B269586BCD299DB2560A0D1359446 21F2B319223D5EEA11EE21AB9254237748D1718365C3306D86FD34756B2FF66F");
    element_set_str(h3, "B5E2403ADE562F5A31ADA4C4521FEEB7897E338A83BA49A354EED1A6355D631 E07BBC42482F395A45276C852C17B7965F6B74B0800D3B6706379B298581E17");

    element_set_str(g1, "9B9F56007ED7AC95D4DBEDD86EFC1615C00BBEB6F0F5BDB0A2292804717FD8F 152B2B9C2DEC1F5E963B3AD6D8E4E5E9E431A2244190711BC38F49359412DB21");
    element_set_str(g2, "1C066BFB76620350F2F68EC29FF49FF51A80A7ABE1B4DA6545DA55D95DBAB46E 172DCA8795031800911AF19DC24107848BDF7CAC7418C611B3CB0A76072DDEB9");
    element_set_str(g3, "D34B41BEA33B724F5E7C194018B1B2EFB28BB248A09D6C985CB5A285023DDAF 14207A55B0C2240D9BF7505F79318A3ABA805AEF1840F72A8F7E50F0A8AEE9DF");

    element_set(((Element *)d1->data)[0], h1);
    element_set(((Element *)d1->data)[1], h2);
    element_set(((Element *)d1->data)[2], h3);

    element_set(((Element *)d2->data)[0], g1);
    element_set(((Element *)d2->data)[1], g2);
    element_set(((Element *)d2->data)[2], g3);

    element_set(((Element *)d->data)[0], d1);
    element_set(((Element *)d->data)[1], d2);

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        bn254_fp12_mul(c, a, b);
    }
    t2 = rdtsc();

    printf("element mul: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //--------------------
    //  sqr
    //--------------------

    element_sqr(c, a);
    element_mul(d, a, a);

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_sqr(c, a);
    }
    t2 = rdtsc();

    printf("element sqr: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //--------------------
    //  random
    //--------------------
    element_random(b);

    //--------------------
    //  inv
    //--------------------
    element_mul(c, a, b);
    element_inv(b, b);
    element_mul(c, c, b);
    element_inv(d, a);
    element_mul(d, a, d);

    assert(element_cmp(c, a) == 0);
    assert(element_is_one(d));

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_inv(c, a);
    }
    t2 = rdtsc();

    printf("element inv: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //--------------------
    //  pow
    //--------------------
    mpz_init(exp);

    mpz_set_str(exp, "103DC9103A8BDD2D08BA8E58EC4AD8AA108177A04CFBBAE6FB7CE07A88E298C3", 16);

    element_pow(c, a, exp);
    element_set_str(d, "1C0084289F0BDD44FAF99FF849E9B8513BBF8A1638E4182B36E47B2DBC03B58A D3E60463F0CDBC861118BE3D448A4F44D0FAF52023351530F56FB9C4E74DEC0 24B9CC7B476D6A9636B2C26117153571CB1510CED2E8A7695BF1372A2BF24E9A 1A62162DE250F132418B7D7193D9AAB496260A4B8C327C6E4A5AF47B6899A1DD 1B2C991268150F299D7F3F66264D8CA5D00DD28C85ED138876CCE3BE6F8683DA FF55E46B7384A0C8EA0F9ABFBFB4CE893EA2E26B6B95B886C766A5C2A40A8AD 56196BA2FB71F737F7F49E515A41AA8A81ACC5DCCD1E7980F107673EA9C1AA7 AC6E66F1AA80D12046178355FA4A7F9C52592857F2F26D0B7A2DAB2ED685E19 BCABC62F9C9B3C8811793FC4DBF136F08E76145ED5244CD8F70B53013A8297A DEA2777027CB2D68EABFE93CE1CC41B9C03C9589E95093AD9E3D689D5081D2D 17CE647195529EF42A7CAAB6F342E2C30D2B1D2F613664FA38502824E9D965BC 1F1D0CF688A1CC8B15BAB38C00C3F8BD7F88BD1A8CA96E73249FE34C15AACF07");

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_pow(b, a, exp);
    }
    t2 = rdtsc();

    printf("element pow with torsion: %.2lf [clock]\n", (double)(t2 - t1) / N);

    mpz_set(exp, f->order);

    for (i = 0; i < 10; i++)
    {
        element_random(a);

        element_pow(b, a, exp);

        assert(element_cmp(b, a) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        element_pow(b, a, exp);
    }
    t2 = rdtsc();

    printf("element pow with order: %.2lf [clock]\n", (double)(t2 - t1) / M);

    mpz_clear(exp);

    //--------------------
    //  clear
    //--------------------
    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(d);

    element_clear(d1);
    element_clear(d2);

    element_clear(h1);
    element_clear(h2);
    element_clear(h3);
    element_clear(g1);
    element_clear(g2);
    element_clear(g3);
}
//============================================
//   test for sqrt
//============================================
void test_sqrt(Field f)
{
    int i;
    unsigned long long int t1, t2;
    Element a, b, c, d;

    element_init(a, f);
    element_init(b, f);
    element_init(c, f);
    element_init(d, f);

    for (i = 0; i < 10; i++)
    {
        element_random(a);
        element_sqr(b, a);

        assert(element_is_sqr(b));

        element_sqrt(c, b);
        element_sqr(d, c);

        assert(element_cmp(d, b) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_is_sqr(b);
    }
    t2 = rdtsc();

    printf("element is sqr: %.2lf [clock]\n", (double)(t2 - t1) / N);

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        element_sqrt(c, b);
    }
    t2 = rdtsc();

    printf("element sqrt: %.2lf [clock]\n", (double)(t2 - t1) / M);

    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(d);
}

//============================================
//   Frobenius Map \phi_p
//============================================
void test_frob(Field f)
{
    int i;
    unsigned long long int t1, t2;
    mpz_t p, p2, p3;
    Element a, b, c;

    mpz_init(p);
    mpz_init(p2);
    mpz_init(p3);

    mpz_set(p, *field_get_char(f));
    mpz_mul(p2, p, p);
    mpz_mul(p3, p2, p);

    element_init(a, f);
    element_init(b, f);
    element_init(c, f);

    for (i = 0; i < 100; i++)
    {
        element_random(a);
        element_pow(b, a, p);
        bn254_fp12_frob_p(c, a);

        assert(element_cmp(b, c) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        bn254_fp12_frob_p(c, a);
    }
    t2 = rdtsc();

    printf("element frob p: %.2lf [clock]\n", (double)(t2 - t1) / N);

    for (i = 0; i < 100; i++)
    {
        element_random(a);
        element_pow(b, a, p2);
        bn254_fp12_frob_p2(c, a);

        assert(element_cmp(b, c) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        bn254_fp12_frob_p2(c, a);
    }
    t2 = rdtsc();

    printf("element frob p2: %.2lf [clock]\n", (double)(t2 - t1) / N);

    for (i = 0; i < 100; i++)
    {
        element_random(a);
        element_pow(b, a, p3);
        bn254_fp12_frob_p3(c, a);

        assert(element_cmp(b, c) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        bn254_fp12_frob_p3(c, a);
    }
    t2 = rdtsc();

    printf("element frob p3: %.2lf [clock]\n", (double)(t2 - t1) / N);

    mpz_clear(p);
    mpz_clear(p2);
    mpz_clear(p3);

    element_clear(a);
    element_clear(b);
    element_clear(c);
}

//============================================
//   i/o test
//============================================
void test_io(Field f)
{
    int i;
    unsigned long long int t1, t2;
    char a_str[780];

    size_t blen;
    unsigned char b_str[384];

    Element a, b, c;

    element_init(a, f);
    element_init(b, f);
    element_init(c, f);

    for (i = 0; i < 100; i++)
    {
        element_random(a);

        element_get_str(a_str, a);
        element_set_str(c, a_str);

        assert(element_cmp(a, c) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_get_str(a_str, a);
    }
    t2 = rdtsc();

    printf("element get string: %.2lf [clock]\n", (double)(t2 - t1) / N);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_set_str(c, a_str);
    }
    t2 = rdtsc();

    printf("element set string: %.2lf [clock]\n", (double)(t2 - t1) / N);

    for (i = 0; i < 100; i++)
    {
        element_random(b);

        element_to_oct(b_str, &blen, b);
        element_from_oct(c, b_str, blen);

        assert(element_cmp(b, c) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_to_oct(b_str, &blen, b);
    }
    t2 = rdtsc();

    printf("element to octet string: %.2lf [clock]\n", (double)(t2 - t1) / N);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_from_oct(c, b_str, blen);
    }
    t2 = rdtsc();

    printf("element from octet string: %.2lf [clock]\n", (double)(t2 - t1) / N);

    element_clear(a);
    element_clear(b);
    element_clear(c);
}

//============================================
//  main program
//============================================
int main(void)
{
    Field fa, fb;

    // test for beuchat's methods
    field_init(fa, "bn254_fp12a");
    test_feature(fa);
    test_arithmetic_operation_beuchat(fa);
    test_sqrt(fa);
    test_frob(fa);
    test_io(fa);

    // test for aranha's methods
    field_init(fb, "bn254_fp12b");
    test_feature(fb);
    test_arithmetic_operation_aranha(fb);
    test_sqrt(fb);
    test_frob(fb);
    test_io(fb);

    field_clear(fa);
    field_clear(fb);

    fprintf(stderr, "ok\n");

    return 0;
}
