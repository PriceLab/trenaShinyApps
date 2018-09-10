import {Selector} from 'testcafe';

fixture `Getting Started`
    .page `http://trena.systemsbiology.net`;

test('My first test', async (t) => {
    await t
        .typeText('#username', 'tesla')
        .typeText('#password', 'password')
        .click('.form-signin > button')
        .expect(Selector("#applist > ul > li:nth-child(7) > a"));
   });

test('My second test', async function foo(t){
    await t
        .typeText('#username', 'tesla')
        .typeText('#password', 'password')
        .click('.form-signin > button')
        .expect(Selector("a").withText("trena INPP5D"));
    });
