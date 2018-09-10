import {Selector} from 'testcafe';

fixture `trena shiny app admin login, check app list`
    .page `http://trena.systemsbiology.net`;

test('login', async (t) => {
    await t
        .typeText('#username', 'paul')
        .typeText('#password', 'leo.paul')
        .click('.form-signin > button')
        //.expect(Selector("#applist > ul > li:nth-child(7) > a"))
        //.expect(Selector("a").withText('trena INPP5D'))
        .expect(Selector("a").withText("cyjShiny demo"));
});

test('FRD3', async t => {
   await t
        .expect(Selector("a").withText("trena FRD3"));
   });


test('test Hello Application', async t => {
    await t
        .expect(Selector("a").withText("Hello Application"));
    });


test('cyjDemo', async t => {
    await t
        .expect(Selector("a").withText("cyjDemo"));
    });

