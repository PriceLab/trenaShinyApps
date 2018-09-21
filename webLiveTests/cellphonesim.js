import {Selector} from 'testcafe';

fixture `Getting Started`
    .page `http://cellphonesim.systemsbiology.net`;

test('launch chiaa app', async function foo(t){

    const pageLoaded = Selector("#phoneTreeButton");
    await t
        .expect(pageLoaded.exists).ok()
    });
